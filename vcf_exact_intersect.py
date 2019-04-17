#!/Users/danielence/anaconda/bin/python
#Daniel Ence January 23, 2019
from cyvcf2 import VCF
import argparse
import sys

def main(args):
    exact_vcf = vcf_exact_intersect(args.a,args.b)
    sys.stderr.write(exact_vcf.get_filenames() + "\n")
    exact_vcf.get_exact_intersection()
    exact_vcf.print_result()

class vcf_exact_intersect(object):

    def __init__(self):
        self.__a_vcf_file=""
        self.__b_vcf_file=""
        self.__intersection_results=[]

    def __init__(self,a_file,b_file):
        self.__a_vcf_file=a_file
        self.__b_vcf_file=b_file
        self.__intersection_results=[]

    def get_filenames(self):
        return "a_file is:\t" + self.__a_vcf_file + "\n" + "b_file is:\t" + self.__b_vcf_file

    def get_exact_intersection(self):
        #make sorted map of positions of VCFs with the variant objects
        a_vcf_map = self.__make_VCF_map(VCF(self.__a_vcf_file))
        b_vcf_map = self.__make_VCF_map(VCF(self.__b_vcf_file))

        sys.stderr.write("this many in \"a\" file: " + str(len(a_vcf_map.keys())) + "\n")
        sys.stderr.write("this many in \"b\" file: " + str(len(b_vcf_map.keys())) + "\n")

        for a_position in a_vcf_map.iterkeys():
            if b_vcf_map.has_key(a_position):
                self.__intersection_results.append(a_vcf_map[a_position])

        sys.stderr.write("this many in the intersection: " + str(len(self.__intersection_results)))

    def __make_VCF_map(self,vcf_obj):

        vcf_map = {}
        for variant in vcf_obj:
            vcf_map.setdefault(self.__make_variant_name(variant),variant)

        return vcf_map

    def __make_variant_name(self,variant_obj):
        return("_".join((variant_obj.CHROM,str(variant_obj.start),str(variant_obj.end))))

    def print_result(self):
        for variant in self.__intersection_results:
            sys.stdout.write(str(variant))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--a",type=str,help="")
    parser.add_argument("--b",type=str,help="")
    args = parser.parse_args()
    main(args)
