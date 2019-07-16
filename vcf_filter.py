#!/usr/bin/python
#Daniel Ence
#

#need to write a filter to get to one SNP per scaffold
#need to write a filter for SNPs that segregate 1:1 in haploid pop or 1:2:1 in diploid pop

import argparse
from cyvcf2 import VCF
from scipy.stats import chisquare
import array from numpy as np.array

class MyVCFFilter(object):
    def __init__(self,vcf_file):
        self.__my_vcf_file = vcf_file
        self.__my_filtered_variants = []
        self.__expected_freqs = []
        #sys.stderr.write("checking expected freqs in constructor\n")
        #sys.stderr.write(str(self.__expected_freqs))
        #sys.stderr.write("\n")

    def filter_segregant(self,alpha=0.05):

        #check to see if we've already started filtering on another criteria
        if(len(self.__my_filtered_variants) == 0):
            #sys.stderr.write("Here once\n")
            curr_VCF_file = VCF(self.__my_vcf_file)
            #for every variant in the vcf file
            for variant in curr_VCF_file:
                #get the non missing genotypes
                #maybe add a filter for missingness in data?
                if(self.properly_segregates(self.get_nonmissing_gts(variant.gt_types), alpha)):
                    self.__my_filtered_variants.append(variant)
                    #sys.stderr.write("This many variants passed:\t" + str(len(self.__my_filtered_variants)))
                    #sys.stderr.write("\n")


        else:
            tmp_variant_list = []
            for variant in self.__my_filtered_variants:
                if (self.properly_segregates(self.get_nonmissing_gts(variant.gt_types), alpha)):
                    tmp_variant_list.append(variant)
            self.__my_filtered_variants = tmp_variant_list

    def set_expected_freqs(self, expected_freqs_string):
        string_freqs = expected_freqs_string.split(",")
        self.__expected_freqs.append(float(string_freqs[0]))
        self.__expected_freqs.append(float(string_freqs[1]))
        #sys.stderr.write("checking expected freqs in setter")
        #sys.stderr.write("\n")
        #sys.stderr.write(str(self.__expected_freqs))
        #sys.stderr.write("\n")

    def get_nonmissing_gts(self, gt_type_list):
        #from brentp's cyvc2 README at: https://github.com/brentp/cyvcf2
        # numpy arrays of specific things we pull from the sample fields.
        #gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT

        #assumming a test for 1:1 segregation in a haploid population for now
        nonmissing_list = []
        for gt in gt_type_list:
            if(gt != 2):
                nonmissing_list.append(gt)
        return nonmissing_list

    def properly_segregates(self, gt_list, alpha):

        #sys.stderr.write("at the top of properly segregates\n")
        #sys.stderr.write(str((gt_list.count(0),gt_list.count(1),gt_list.count(3))))
        #sys.stderr.write("\n")

        if(gt_list.count(1) == 0):
            #if there are no hets, then test 1:1
            pvalue = self.__do_one_to_one_chi2_test(gt_list)
        elif(gt_list.count(3) == 0):
            pvalue = self.__do_one_to_one_chi2_test(gt_list)
        elif(gt_list.count(0) == 0):
            pvalue = self.__do_one_to_one_chi2_test(gt_list)
        else:
            #assuming Diploid het. cross here
            pvalue = self.__do_one_two_one_chi2_test(gt_list)

        #sys.stderr.write("Checking the alpha and pvalue variablees\n")
        #sys.stderr.write(str(alpha))
        #sys.stderr.write("\n")
        #sys.stderr.write(str(pvalue[1]))
        #sys.stderr.write("\n")
        #sys.stderr.write(str(pvalue[1] > alpha))
        #sys.stderr.write("\n")

        return pvalue[1] > alpha

    def __do_one_to_one_chi2_test(self, gt_list):
        #sys.stderr.write("one to one test")

        if(len(self.__expected_freqs) > 0):

            expected_counts = []
            for expected_freq in self.__expected_freqs:
                expected_counts.append(round(expected_freq * len(gt_list)))
            #sys.stderr.write("****************Expected counts**************\n")
            #sys.stderr.write(str(expected_counts))
            #sys.stderr.write("\n")

            if(gt_list.count(1) == 0):
                return chisquare(f_obs=[gt_list.count(0), gt_list.count(3)], f_exp=expected_counts)
            elif(gt_list.count(3) == 0):
                return chisquare(f_obs=[gt_list.count(0), gt_list.count(1)], f_exp=expected_counts)
            elif(gt_list.count(0) == 0):
                return chisquare(f_obs=[gt_list.count(1), gt_list.count(3)], f_exp=expected_counts)
        else:
            if (gt_list.count(1) == 0):
                return chisquare(f_obs=[gt_list.count(0), gt_list.count(3)])
            elif (gt_list.count(3) == 0):
                return chisquare(f_obs=[gt_list.count(0), gt_list.count(1)])
            elif (gt_list.count(0) == 0):
                return chisquare(f_obs=[gt_list.count(1), gt_list.count(3)])

    def __do_one_two_one_chi2_test(self, gt_list):

        #sys.stderr.write("one to two to one test")
        #sys.stderr.write(str((gt_list.count(0), gt_list.count(1), gt_list.count(3))))
        if (len(self.__expected_freqs) > 0):

            expected_counts = []
            for expected_freq in self.__expected_freqs:
                expected_counts.append(round(expected_freq * len(gt_list)))
            #sys.stderr.write(str(expected_counts))

            return chisquare(f_obs=[gt_list.count(0) + (gt_list.count(1) / 2),
                                    gt_list.count(3) + (gt_list.count(1) / 2)], f_exp=expected_counts)
        else:
            return chisquare(f_obs=[gt_list.count(0) + (gt_list.count(1) / 2),
                                    gt_list.count(3) + (gt_list.count(1) / 2)])

    def dump_header(self):
        for line in self.__header_lines:
            print(line)

    def dump_filtered_vars(self):
        #need to figure out how to print the header
        #print(str(var) for var in self.__my_filtered_variants)
        #sys.stderr.write(str(len(self.__my_filtered_variants)))
        #sys.stderr.write("\n")
        for var in self.__my_filtered_variants:
            print(str(var).strip())

    def filter_missing(self,missing_max=0.2):
        if(len(self.__my_filtered_variants) == 0):
            curr_VCF_file = VCF(self.__my_vcf_file)

            for variant in curr_VCF_file:
                if(self.__passes_missing(variant.gt_types, missing_max)):
                    self.__my_filtered_variants.append(variant)
        else:
            tmp_variant_list = []
            for variant in self.__my_filtered_variants:
                if(self.__passes_missing(variant.gt_types, missing_max)):
                    tmp_variant_list.append(variant)
            self.__my_filtered_variants = tmp_variant_list

    def __passes_missing(self, gt_list, missing_max):
        return self.__calc_missingness(gt_list) <= missing_max


    def filter_SNPs_per_scaffold(self, num_SNPs_per_scaffold, method):

        if(len(self.__my_filtered_variants) == 0):
            curr_VCF_file = VCF(self.__my_vcf_file)

            tmp_var_list = []
            for var in curr_VCF_file:
                tmp_var_list.append(var)
            self.__filter_SNPs_per_scaffold(tmp_var_list, num_SNPs_per_scaffold, method)
        else:
            self.__filter_SNPs_per_scaffold(self.__my_filtered_variants, num_SNPs_per_scaffold, method)

    def __filter_SNPs_per_scaffold(self, var_list, num_per_scaffold, method):

        #get variants that share scaffolds. Assumes sorted by scaffold in the vcf file.
        curr_scaff_variants = []
        chosen_variants = []
        for variant in var_list:
            #if it is the first one, just add it to the list
            if(len(curr_scaff_variants) == 0):
                curr_scaff_variants.append(variant)
            #if it isnot the first one, check if it matches the CHROM of the list
            elif(variant.CHROM == curr_scaff_variants[0].CHROM):
                #if yes, just add it to the list
                curr_scaff_variants.append(variant)
            else:
                #if not, then we have a group of SNPs to filter
                chosen_variants.extend(self.__choose_SNPs(curr_scaff_variants, num_per_scaffold, method))
                curr_scaff_variants = [variant]
        self.__my_filtered_variants = chosen_variants

    def __choose_SNPs(self, list_of_variants, num_to_choose, method):

        #return list_of_variants[0]
        if(method == "missing"):
            return self.__choose_N_SNPs_missingness(list_of_variants, num_to_choose)
        #elif(method == "max_freq"):
        #    return self.__choose_N_SNPs_max_freq(list_of_variants, num_to_choose)
        #elif(method == "min_freq"):
        #    return self.__choose_N_SNPs_min_freq(list_of_variants, num_to_choose)
        else:
        #    #default is to choose ones closest to the front of the scaffold
            return self.__choose_N_SNPs_min_scaff_pos(list_of_variants, num_to_choose)

    def __choose_N_SNPs_missingness(self, list_of_variants, num_to_choose):

        var_to_missing_map = {}
        for var in list_of_variants:
            missingness = self.__calc_missingness(var.gt_types)
            var_to_missing_map.setdefault(var, missingness)

        import operator
        sorted_map = sorted(var_to_missing_map.items(), key=operator.itemgetter(1))

        chosen_list = []
        for i in range(num_to_choose):
            chosen_list.append(sorted_map[i][0])

        return chosen_list

    def __calc_missingness(self, gt_list):

        missing_count = 0
        for gt in gt_list:
            if(gt == 2):
                missing_count = missing_count + 1
        return float(missing_count) / float(len(gt_list))

    def haplo_diplo_filter(self):
        vcf = VCF(self.__my_vcf_file)

        for variant in vcf:
            #list comprehension. Neato!
            print(variant.genotypes)
            #tmp_genotypes = np.array([2 if x==1 else x for x in variant.genotypes])
            tmp_genotypes = [2 if x==1 else x for x in variant.genotypes]
            variant.genotypes = tmp_genotypes
            self.__my_filtered_variants.append(variant)

def main(args):
    filter_obj = MyVCFFilter(args.vcf_file)
    #first filter on missing data
    if(args.haplo_diplo_missing):
        filter_obj.haplo_diplo_filter()

    if(args.filter_missing==True):
        if(args.missing_max):
            filter_obj.filter_missing(args.missing_max)
        else:
            filter_obj.filter_missing()

    #then filter on segregation
    if(args.filter_chi==True):
        if(args.expected_freqs):
            filter_obj.set_expected_freqs(args.expected_freqs)

        if(args.chi_square_alpha):
                filter_obj.filter_segregant(args.chi_square_alpha)
        else:
            filter_obj.filter_segregant()

    if(args.number_of_SNPs_per_scaffold):
        if(args.SNP_choose_method):
            filter_obj.filter_SNPs_per_scaffold(args.number_of_SNPs_per_scaffold,args.SNP_choose_method)
        else:
            filter_obj.filter_SNPs_per_scaffold(args.number_of_SNPs_per_scaffold)

    #need to figure out how to dump the header info
    #dump filtered vcf
    filter_obj.dump_filtered_vars()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file",type=str,help="vcf file that we're going to filter")
    parser.add_argument("--haplo_diplo_missing",help="tells the script this is haploid called as diploid file. change hets to missing",action='store_true')
    parser.add_argument("--filter_missing",help="tells the script to filter on missing data",action='store_true')
    parser.add_argument("--missing_max",type=float,help="maximum percentage of missing data")
    parser.add_argument("--filter_chi",help="tells the script to filter with a chisquare test for proper segregation",action='store_true')
    parser.add_argument("--expected_freqs",help="a string delimited list of floats to specify your own expected frequencies for the chi2 teset")
    parser.add_argument("--chi_square_alpha",type=float,help="alpha level for the chisquare test")
    #maybe add an option to test for different segregation types?
    #add an option to filter on missingness
    parser.add_argument("--number_of_SNPs_per_scaffold",type=int,help="number of SNPs")
    parser.add_argument("--SNP_choose_method",default="missing", type=str,help="method to choose SNPs per scaffold, one of \'missing\', \'max_freq\',\'min_freq\', default is first on the scaffold")

    args = parser.parse_args()
    main(args)
