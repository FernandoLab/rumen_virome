#! /usr/bin/env python

"""Split deep viral samples from SRA into their runs for better QC
processing."""

from Bio import SeqIO
import itertools
from multiprocessing import Process


def split_v4_r1():
    v4_r1_iterator = SeqIO.parse(open('raw_deep_viral/V4_illumina_1.fastq', 'rU'),
            "fastq")

    run2_iterator = itertools.islice(v4_r1_iterator, 557258)
    output_handle = open('raw_deep_viral/V4_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run1_iterator = itertools.islice(v4_r1_iterator, 1235638)
    output_handle = open('raw_deep_viral/V4_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = v4_r1_iterator
    output_handle = open('raw_deep_viral/V4_illumina_run4_R1.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v4_r2():
    v4_r2_iterator = SeqIO.parse(open('raw_deep_viral/V4_illumina_2.fastq', 'rU'),
            "fastq")

    run2_iterator = itertools.islice(v4_r2_iterator, 557258)
    output_handle = open('raw_deep_viral/V4_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run1_iterator = itertools.islice(v4_r2_iterator, 1235638)
    output_handle = open('raw_deep_viral/V4_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = v4_r2_iterator
    output_handle = open('raw_deep_viral/V4_illumina_run4_R2.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v9_r1():
    v9_r1_iterator = SeqIO.parse(open('raw_deep_viral/V9_illumina_1.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v9_r1_iterator, 123076)
    output_handle = open('raw_deep_viral/V9_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v9_r1_iterator, 1769872)
    output_handle = open('raw_deep_viral/V9_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v9_r1_iterator
    output_handle = open('raw_deep_viral/V9_illumina_run3_R1.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v9_r2():
    v9_r2_iterator = SeqIO.parse(open('raw_deep_viral/V9_illumina_2.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v9_r2_iterator, 123076)
    output_handle = open('raw_deep_viral/V9_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v9_r2_iterator, 1769872)
    output_handle = open('raw_deep_viral/V9_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v9_r2_iterator
    output_handle = open('raw_deep_viral/V9_illumina_run3_R2.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v10_r1():
    v10_r1_iterator = SeqIO.parse(open('raw_deep_viral/V10_illumina_1.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v10_r1_iterator, 36861)
    output_handle = open('raw_deep_viral/V10_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = itertools.islice(v10_r1_iterator, 2552876)
    output_handle = open('raw_deep_viral/V10_illumina_run4_R1.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v10_r1_iterator, 3553281)
    output_handle = open('raw_deep_viral/V10_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v10_r1_iterator
    output_handle = open('raw_deep_viral/V10_illumina_run3_R1.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v10_r2():
    v10_r2_iterator = SeqIO.parse(open('raw_deep_viral/V10_illumina_2.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v10_r2_iterator, 36861)
    output_handle = open('raw_deep_viral/V10_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = itertools.islice(v10_r2_iterator, 2552876)
    output_handle = open('raw_deep_viral/V10_illumina_run4_R2.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v10_r2_iterator, 3553281)
    output_handle = open('raw_deep_viral/V10_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v10_r2_iterator
    output_handle = open('raw_deep_viral/V10_illumina_run3_R2.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v11_r1():
    v11_r1_iterator = SeqIO.parse(open('raw_deep_viral/V11_illumina_1.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v11_r1_iterator, 261550)
    output_handle = open('raw_deep_viral/V11_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v11_r1_iterator, 5152548)
    output_handle = open('raw_deep_viral/V11_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = v11_r1_iterator
    output_handle = open('raw_deep_viral/V11_illumina_run4_R1.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v11_r2():
    v11_r2_iterator = SeqIO.parse(open('raw_deep_viral/V11_illumina_2.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v11_r2_iterator, 261550)
    output_handle = open('raw_deep_viral/V11_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v11_r2_iterator, 5152548)
    output_handle = open('raw_deep_viral/V11_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = v11_r2_iterator
    output_handle = open('raw_deep_viral/V11_illumina_run4_R2.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v12_r1():
    v12_r1_iterator = SeqIO.parse(open('raw_deep_viral/V12_illumina_1.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v12_r1_iterator, 76336)
    output_handle = open('raw_deep_viral/V12_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v12_r1_iterator, 2617624)
    output_handle = open('raw_deep_viral/V12_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = itertools.islice(v12_r1_iterator, 1846464)
    output_handle = open('raw_deep_viral/V12_illumina_run4_R1.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v12_r1_iterator
    output_handle = open('raw_deep_viral/V12_illumina_run3_R1.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v12_r2():
    v12_r2_iterator = SeqIO.parse(open('raw_deep_viral/V12_illumina_2.fastq', 'rU'),
            "fastq")

    run1_iterator = itertools.islice(v12_r2_iterator, 76336)
    output_handle = open('raw_deep_viral/V12_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v12_r2_iterator, 2617624)
    output_handle = open('raw_deep_viral/V12_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run4_iterator = itertools.islice(v12_r2_iterator, 1846464)
    output_handle = open('raw_deep_viral/V12_illumina_run4_R2.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v12_r2_iterator
    output_handle = open('raw_deep_viral/V12_illumina_run3_R2.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v13_r1():
    v13_r1_iterator = SeqIO.parse(open('raw_deep_viral/V13_illumina_1.fastq', 'rU'),
            "fastq")

    run4_iterator = itertools.islice(v13_r1_iterator, 442696)
    output_handle = open('raw_deep_viral/V13_illumina_run4_R1.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v13_r1_iterator, 195176)
    output_handle = open('raw_deep_viral/V13_illumina_run2_R1.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run1_iterator = itertools.islice(v13_r1_iterator, 210184)
    output_handle = open('raw_deep_viral/V13_illumina_run1_R1.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v13_r1_iterator
    output_handle = open('raw_deep_viral/V13_illumina_run3_R1.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def split_v13_r2():
    v13_r2_iterator = SeqIO.parse(open('raw_deep_viral/V13_illumina_2.fastq', 'rU'),
            "fastq")

    run4_iterator = itertools.islice(v13_r2_iterator, 442696)
    output_handle = open('raw_deep_viral/V13_illumina_run4_R2.fastq', 'w')
    SeqIO.write(run4_iterator, output_handle, 'fastq')
    output_handle.close()

    run2_iterator = itertools.islice(v13_r2_iterator, 195176)
    output_handle = open('raw_deep_viral/V13_illumina_run2_R2.fastq', 'w')
    SeqIO.write(run2_iterator, output_handle, 'fastq')
    output_handle.close()

    run1_iterator = itertools.islice(v13_r2_iterator, 210184)
    output_handle = open('raw_deep_viral/V13_illumina_run1_R2.fastq', 'w')
    SeqIO.write(run1_iterator, output_handle, 'fastq')
    output_handle.close()

    run3_iterator = v13_r2_iterator
    output_handle = open('raw_deep_viral/V13_illumina_run3_R2.fastq', 'w')
    SeqIO.write(run3_iterator, output_handle, 'fastq')
    output_handle.close()

def run_parallel(*fns):
  proc = []
  for fn in fns:
    p = Process(target=fn)
    p.start()
    proc.append(p)
  for p in proc:
    p.join()

if __name__ == '__main__':
    run_parallel(split_v4_r1, split_v4_r2, split_v9_r1, split_v9_r2,
            split_v10_r1, split_v10_r2, split_v11_r1, split_v11_r2,
            split_v12_r1, split_v12_r2, split_v13_r1, split_v13_r2)
