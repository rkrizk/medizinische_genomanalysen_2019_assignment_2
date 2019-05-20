#!/usr/bin/python3
import argparse
import sys
import vcf

parser = argparse.ArgumentParser(description='Summary of VCF file contents')

parser.add_argument('-1', type=str, dest='file_1',
                    help='path to first VCF file')
parser.add_argument('-2', type=str, dest='file_2',
                    help='path to second VCF file')


if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
file_1 = args.file_1
file_2 = args.file_2

__author__ = 'Richard Kriz'


class Assignment2:
    
    def __init__(self):
        print("PyVCF version: %s" % vcf.VERSION)

    def get_average_quality_of_file(self):
        """
        Get the average PHRED quality of all variants
        :return:
        """
        var_counter = 0
        quality_sum = 0
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                quality_sum += record.samples[0]['GQ']
                var_counter += 1
        return round(quality_sum / var_counter, 2)

    def get_total_number_of_variants_of_file(self):
        """
        Get the total number of variants
        :return: total number of variants
        """
        var_counter = 0
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                var_counter += 1
        return var_counter

    def get_variant_caller_of_vcf(self):
        """
        Return the variant caller name
        :return: 
        """
        var_caller_list = []
        with open(file_1) as my_vcf_fh:
            for line in my_vcf_fh:
                if not line.startswith("#"):
                    query = line[line.find("callsetnames")+13:]
                    query = query[:query.find(";")]
                    query = query.split(",")
                    for var_caller in query:
                        if var_caller not in var_caller_list:
                            var_caller_list.append(var_caller)
        return var_caller_list
        
    def get_human_reference_version(self):
        """
        Return the genome reference version
        :return: 
        """
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            try:
                reference = vcf_reader.metadata['reference']
            except:
                reference = "Not specified"

        return reference

    def get_number_of_indels(self):
        """
        Return the number of identified INDELs
        :return:
        """
        indel_counter = 0
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                if record.is_indel:
                    indel_counter += 1
        return indel_counter

    def get_number_of_snvs(self):
        """
        Return the number of SNVs
        :return: 
        """
        snv_counter = 0
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                if record.is_snp:
                    snv_counter += 1
        return snv_counter
        
    def get_number_of_heterozygous_variants(self):
        """
        Return the number of heterozygous variants
        :return: 
        """
        het_counter = 0
        heterozygous = ["1|0", "0|1", "1/0", "0/1"]
        with open(file_1) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for record in vcf_reader:
                if record.samples[0]['GT'] in heterozygous:
                    het_counter += 1
        return het_counter

    def merge_chrs_into_one_vcf(self):
        """
        Creates one VCF containing all variants of chr21 and chr22
        :return:
        """

        reader_1 = vcf.Reader(open(file_1))
        writer_1 = vcf.Writer(open("merge.vcf", "w+"), reader_1)
        reader_2 = vcf.Reader(open(file_1))
        writer_2 = vcf.Writer(open("merge.vcf", "a"), reader_2)
        for record in reader_1:
            writer_1.write_record(record)

        for record in reader_2:
            writer_2.write_record(record)

        line_counter = 0
        with open("merge.vcf") as merge_file:
            for line in merge_file:
                line_counter += 1

        return line_counter

    def print_summary(self):
        print("\nAverage Quality of file: {}".format(self.get_average_quality_of_file()))
        print("Total number of variants in file: {}".format(self.get_total_number_of_variants_of_file()))
        print("Variant Caller: {}".format(self.get_variant_caller_of_vcf()))
        print("Reference Version: {}".format(self.get_human_reference_version()))
        print("Number of INDELs: {}".format(self.get_number_of_indels()))
        print("Number of SNVs: {}".format(self.get_number_of_snvs()))
        print("Number of heterozygous variants: {}".format(self.get_number_of_heterozygous_variants()))
        print("Number of variants in merged file: {}\n".format(self.merge_chrs_into_one_vcf()))


def main():
    print("Assignment 2")
    assignment2 = Assignment2()
    assignment2.print_summary()
    print("Done with assignment 2")
        
        
if __name__ == '__main__':
    main()
