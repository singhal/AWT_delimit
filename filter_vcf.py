import re
import gzip

# raw_vcf = '/scratch/drabosky_flux/sosi/AWT_delimit/reference_bias/Lampropholis_coggeri_S.raw2.vcf'
raw_vcf = '/scratch/drabosky_flux/sosi/AWT_delimit/original_approach/variants/Lampropholis_coggeri_S.qual_filtered.cov_filtered.vcf.gz'
# filt_vcf = '/scratch/drabosky_flux/sosi/AWT_delimit/variants/Lampropholis_coggeri_S.qual_filtered.cov_filtered.qual.vcf.gz'
filt_vcf = ''

def run(raw_vcf, filt_vcf):
        f = gzip.open(raw_vcf, 'r')
        # o = gzip.open(filt_vcf, 'w')

        for l in f:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
                        # check if indel
                        # alleles = [d[3]] + re.split(',', d[4])
                        # indel = False
                        # for a in alleles:
                        #        if len(a) > 1:
                        #                indel = True
                        
                        # depth cover
                        # n_inds = len(d[9:])
                        # dp = 10 * n_inds
                        # snp_depth = int(re.search('DP=(\d+)', d[7]).group(1))
                        
			qd = 10
			if re.search('QD=([\d|.]+)', d[7]):
				qd = float(re.search('QD=([\d|.]+)', d[7]).group(1))
			if qd >= 5:
				print(l.rstrip())

run(raw_vcf, filt_vcf)
