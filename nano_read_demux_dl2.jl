# Author: Konrad Herbst, Dan Lou 
# Heidelberg University, ZMBH, Knop Lab
# 2025-09-02
# In-house script for demultiplexing and trimming


using ArgParse
using BioSequences    ### the version has to be v2.0.1
using FASTX
using CodecZlib
using Statistics
using CSV
using DataFrames      

## PARAMETERS
ACCEPTED_ERRORS = 3                                  # maximal numbers of errors per 10 bp allowed for constant_seq search
TN5_ADAPTER = dna"TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" # DNA sequence of the Tn5 adapter used
BC_ERRORS = 1                                        # maximal numbers of allowed errors in BC demux; TODO generalize with ACCEPTED_ERRORS
BC_DISTANCEMETRIC = :levenshtein                     # distance metric used for BC demultiplexing

## helper functions
function newPath(path::AbstractString, suffix::AbstractString, dir = nothing)
    (path, ext) = splitext(path)
    (pathdir, file) = splitdir(path)
    dir = dir == nothing ? pathdir : dir  #the statement sets dir to pathdir if it was initially nothing, otherwise it leaves it unchanged.
    if !isdir(dir)
        # throw(ErrorException("OUTDIR does not exist."))
        println(string("'", dir, "' does not exist. Creating."))
        mkpath(dir)
    end
    return string(dir, "/", file, suffix, ext)
end

function splitextex(path::AbstractString)
    next = ""
    (pathdir, file) = splitdir(path)
    (file, ext) = splitext(file)
    while ext != ""
        next = string(ext, next)
        (file, ext) = splitext(file)
    end
    return (file, next)
end

function format_count(x::Int64, total::Int64)
    return "$x\t($(round(x/total*1000)/10) %)"
end

function parse_commandline()
    s = ArgParseSettings("Trimming and demux script for Tn5-Anchor-Seq protocol (MinION)\n" *
                         "Konrad Herbst, k.herbst@zmbh.uni-heidelberg.de",
                         version = "Version 0.2",
                         add_version = true)

    @add_arg_table! s begin
        "--var-extract", "-x"
            help = "extract given number of nt after constant_seq, write into readname and strip off from R1 and R2; 0=OFF"
            arg_type = Int
            default = 0
        "--minimal-length", "-l"
            help  = "minimal length for the trimmed reads to be retained"
            arg_type = Int
            default = 20
        "--compress", "-z"
            help  = "compress output files. Always done for compressed input files."
            action = :store_true
        "reads"
            help = "Input fastq file."
            required = true
        "samples"
            help = "Samplesheet (CSV format) (with demultiplexing) _OR_ constant_seq (without demultiplexing). " *
                   "Will be decided, based on if samplesheet is a proper file."
            required = true
        "outdir"
            help = "directory to write trimmed reads into"
            required = true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()


stats = Dict("TOTAL" => 0,
             "TRIMMED" => 0,
             "NO_CONSTANT" => 0,
             "NO_TN5ADAPTER" => 0,
             "TOO_SHORT" => 0,
             "NOT_ASSIGNED" => 0)

## determine if sample-description was given via constant_seq (no demux) or samplesheet
demux = isfile(parsed_args["samples"])
if !demux
    constant_seq = parsed_args["samples"]
else
    ## transform sample description file into useful data structures
    samples = CSV.read(parsed_args["samples"], DataFrame)  ###2
    constant_seq = unique(samples[!, :constant_sequence])
    if length(constant_seq) > 1
        throw(ErrorException("Handling of more than one 'constant_seq' not yet implemented."))
    else
        constant_seq = constant_seq[1]
    end
    BC_length = unique([length(x) for x in samples[!, :barcodeRV_sequence]])
    if length(BC_length) > 1
        throw(ErrorException("Differing barcode length not supported."))
    else
        BC_length = BC_length[1]
    end
    ## construct demuxer
    BCs_read = [reverse_complement(LongDNASeq(x)) for x in samples[!, :barcodeRV_sequence]]
    dplxr_read = Demultiplexer(BCs_read, n_max_errors = BC_ERRORS, distance = BC_DISTANCEMETRIC)
end
constant_seq_pattern = ApproximateSearchQuery(LongDNASeq(constant_seq))
adapter_pattern = ApproximateSearchQuery(reverse_complement(TN5_ADAPTER))

## initiate output file handles
(reads_name, reads_ext) = splitextex(parsed_args["reads"])
if !demux
    if splitext(reads_ext)[2] == ".gz" || parsed_args["compress"]
        fhs = FASTQ.Writer(GzipCompressorStream(open(newPath(reads_name * ".fastq.gz", "", parsed_args["outdir"]), "w")))
    else
        fhs = open(FASTQ.Writer, newPath(reads_name * ".fastq", "", parsed_args["outdir"]))
    end
else
    fhs = Vector{FASTQ.Writer}(undef, length(BCs_read))
    for ib in 1:length(BCs_read)
        fname = samples[ib, :sample]           # ":sample" means the column "sample" in samples dataframe
        if splitext(reads_ext)[2] == ".gz" || parsed_args["compress"]
            fhs[ib] = FASTQ.Writer(GzipCompressorStream(open(newPath(fname * ".fastq.gz", "", parsed_args["outdir"]), "w")))
        else
            fhs[ib] = open(FASTQ.Writer, newPath(fname * ".fastq", "", parsed_args["outdir"]))
        end
    end
end
## in addition create a io sink for undetermined sequences
if splitext(reads_ext)[2] == ".gz" || parsed_args["compress"]
    fhs_undet = FASTQ.Writer(GzipCompressorStream(open(newPath("undetermined.fastq.gz", "", parsed_args["outdir"]), "w")))
else
    fhs_undet = open(FASTQ.Writer, newPath("undetermined.fastq", "", parsed_args["outdir"]))
end

## read trimming function which does all the work
function trim_read(read::FASTQ.Record, constant_seq_pattern::ApproximateSearchQuery{LongSequence{DNAAlphabet{4}}},
                   var_extract::Int64 = 0, minlengths::Int64 = 20)
    ## collect fail reasons and initialize demux index
    fail = Nothing
    bci = -1

    seq = FASTQ.sequence(read)
    qua = FASTQ.quality(read)
    ide = FASTQ.identifier(read)
    des = FASTQ.description(read)

    (constant, variable, bc) = ("", "", "")

    if (length(seq) >= (length(constant_seq_pattern.seq) + minlengths))
        # search for end of constant_seq in seq
        ce = last(approxrsearch(seq, constant_seq_pattern, Int(ACCEPTED_ERRORS * round(length(constant_seq_pattern.seq)/10))))
        if ce != -1 # read contains end of constant_seq
            ## search for end of TN5_ADAPTER in seq
            ae = last(approxrsearch(seq, adapter_pattern, Int(ACCEPTED_ERRORS * round(length(adapter_pattern.seq)/10))))
            if ae != -1 # read contains end of TN5_ADAPTER
                ## search for start of TN5_ADAPTER in seq
                as = first(approxsearch(seq, adapter_pattern, Int(ACCEPTED_ERRORS * round(length(adapter_pattern.seq)/10))))
                constant = seq[1:ce]
                seq = seq[(ce+1):(as-1)]
                qua = qua[(ce+1):(as-1)]
                    ## extract variable part
                if var_extract > 0
                    variable = seq[1:min(var_extract, length(seq))]
                    seq = !isempty(seq) ? seq[min(var_extract+1, length(seq)+1):end] : seq
                        qua = !isempty(qua) ? qua[min(var_extract+1, length(qua)+1):end] : qua
                    ide = string(ide, ':', variable)
                end
                des = string(des, ':', constant, ':', bc)
            else
                fail = "NO_TN5ADAPTER"
            end # ae != -1
        else
            fail = "NO_CONSTANT"
        end # ce != -1
    else
        fail = "TOO_SHORT"
    end # length(seq)

    return (FASTQ.Record(ide, des, seq, qua), fail)
end ## trim_demux_read

## read trimming+demux function which does all the work
function trim_demux_read(read::FASTQ.Record, constant_seq_pattern::ApproximateSearchQuery{LongSequence{DNAAlphabet{4}}},
                   dplxr::Demultiplexer{LongSequence{DNAAlphabet{4}}}, var_extract::Int64 = 0, minlengths::Int64 = 20)
    ## collect fail reasons and initialize demux index
    fail = Nothing
    bci = -1

    seq = FASTQ.sequence(read)
    qua = FASTQ.quality(read)
    ide = FASTQ.identifier(read)
    des = FASTQ.description(read)

    (constant, variable, bc) = ("", "", "")

    if (length(seq) >= (length(constant_seq_pattern.seq) + minlengths))
        # search for end of constant_seq in seq
        ce = last(approxrsearch(seq, constant_seq_pattern, Int(ACCEPTED_ERRORS * round(length(constant_seq_pattern.seq)/10))))
        if ce != -1 # read contains end of constant_seq
            ## search for end of TN5_ADAPTER in seq
            ae = last(approxrsearch(seq, adapter_pattern, Int(ACCEPTED_ERRORS * round(length(adapter_pattern.seq)/10)))) 
            if ae != -1 # read contains end of TN5_ADAPTER
                ## try to demultiplex
                if (length(seq) >= (ae + BC_length))
                    bc = seq[(ae+1):(ae+BC_length)]
                    bci = first(demultiplex(dplxr_read, bc))
                    if bci > 0 # read could be demultiplexed
                        ## search for start of TN5_ADAPTER in seq
                        as = first(approxsearch(seq, adapter_pattern, Int(ACCEPTED_ERRORS * round(length(adapter_pattern.seq)/10))))
                        constant = seq[1:ce]
                        seq = seq[(ce+1):(as-1)]
                        qua = qua[(ce+1):(as-1)]
                        ## extract variable part
                        if var_extract > 0
                            variable = seq[1:min(var_extract, length(seq))]
                            seq = !isempty(seq) ? seq[min(var_extract+1, length(seq)+1):end] : seq
                            qua = !isempty(qua) ? qua[min(var_extract+1, length(qua)+1):end] : qua
                            ide = string(ide, ':', variable)
                        end
                        des = string(des, ':', constant, ':', bc)
                    else
                        fail = "NOT_ASSIGNED"
                    end # bci != -1
                else
                    fail = "TOO_SHORT"
                end
            else
                fail = "NO_TN5ADAPTER"
            end # ae != -1
        else
            fail = "NO_CONSTANT"
        end # ce != -1
    else
        fail = "TOO_SHORT"
    end # length(seq)

    return (FASTQ.Record(ide, des, seq, qua), bci, fail)
end

## iterate over handles process reads
if splitext(parsed_args["reads"])[2] == ".gz"
    reads = FASTQ.Reader(GzipDecompressorStream(open(parsed_args["reads"], "r")))
else
    reads = open(FASTQ.Reader, parsed_args["reads"])
end

## process reads
while !eof(reads)
    r = FASTQ.Record()
    read!(reads, r)
    stats["TOTAL"] += 1

    if !demux
        (tr, fail) = trim_read(r, constant_seq_pattern, parsed_args["var-extract"], parsed_args["minimal-length"])
    else
        (tr, samplei, fail) = trim_demux_read(r, constant_seq_pattern, dplxr_read, parsed_args["var-extract"], parsed_args["minimal-length"])
    end

    ## try reverse complement if failed
    if fail == "NO_CONSTANT"
        r = FASTQ.Record(FASTQ.identifier(r), FASTQ.description(r), reverse_complement(FASTQ.sequence(r)), FASTQ.quality(r)[end:-1:1])
        if !demux
            (tr, fail) = trim_read(r, constant_seq_pattern, parsed_args["var-extract"], parsed_args["minimal-length"])
        else
            (tr, samplei, fail) = trim_demux_read(r, constant_seq_pattern, dplxr_read, parsed_args["var-extract"], parsed_args["minimal-length"])
        end
    end # fail

    ## write reads if they are still long enough
    if fail == Nothing
        if length(FASTQ.sequence(tr)) < parsed_args["minimal-length"]
            stats["TOO_SHORT"] += 1
        else
            if !demux
                write(fhs, tr)
            else
                write(fhs[samplei], tr)
            end
            stats["TRIMMED"] += 1
        end ## if tr
    else
        stats[fail] += 1
        ## write to undetermined reads
        r = FASTQ.Record(FASTQ.identifier(r), string(FASTQ.description(r), ':', fail), FASTQ.sequence(r), FASTQ.quality(r))
        write(fhs_undet, tr)
    end # fail

    ## report progress once in a while
    if (stats["TOTAL"] > 0) & (stats["TOTAL"] % 100000 == 0)
        println(stderr, string(stats["TOTAL"], " reads processed."))
    end
end

## report final statistics
println(stderr, string("Results are in: ", parsed_args["outdir"]))
println(stderr, string(stats["TOTAL"], "\treads were read."))
println(stderr, string(format_count(stats["NO_CONSTANT"], stats["TOTAL"]), " reads in which 'constant_seq' could not be identified."))
println(stderr, string(format_count(stats["NO_TN5ADAPTER"], stats["TOTAL"]), " reads in which 'TN5_ADAPTER' could not be identified."))
if demux
    println(stderr, string(format_count(stats["NOT_ASSIGNED"], stats["TOTAL"]), " reads which could not be assigned based on 'samples'."))
end
println(stderr, string(format_count(stats["TOO_SHORT"], stats["TOTAL"]), " reads which were too short after trimming."))
println(stderr, string(format_count(stats["TRIMMED"], stats["TOTAL"]), " reads were trimmed in total."))

## destruct file handles
close(reads)
if !demux
    close(fhs)
else
    for fhss in fhs
        close(fhss)
    end
end
close(fhs_undet)
