rule fastqc:
    input:
        "A.fastq"
    output:
        html="A.html",
        zip="A_fastqc.zip"
    params: "--quiet"
    log:
        "A.log"
    threads: 2
    script:
        "fastqc.py"
