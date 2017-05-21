import sys
import pysam

chromToUse = sys.argv[1]  # 0 for all chromosomes
norm_hetpsns = sys.argv[2]
bam_file = sys.argv[3]
ref_file = sys.argv[4]
base_quality = int(sys.argv[5])
map_quality = int(sys.argv[6])
vcf_quality = int(sys.argv[7])
positions = {}

#  add (position,depth) from the normal hetpositions input file to a dictionary of lists
#    indexed by chromosome
for line in open(norm_hetpsns):
    if not line.strip().startswith("#"):
        chrom = line.split()[0] 
        if chrom == chromToUse or chromToUse == 0:
            position = int(line.strip().split()[1])
            ref_base = line.strip().split()[3]
            nref_base = line.strip().split()[4]
            qual = line.strip().split()[5]
            depth = line.split()[7].split(';')[0].replace('DP=', '')
            genoKey = line.split()[8].split(':')
            geno = line.split()[9].split(':')
            position_data = position, depth, ref_base, nref_base, qual, genoKey, geno
            if chrom not in positions:
                positions[chrom] = []
            positions[chrom].append(position_data)

sample = pysam.AlignmentFile(bam_file)
#reference = pysam.FastaFile(ref_file)
## print header ##
print "Chr\tPosition\tRef\tRefCount\tNref\tNrefCount\tNormQuality\tNormGenotype\tNormPhaseSet"

for chrom in positions:
    i = 0
    for position_data in positions[chrom]:
			position = int(position_data[0])
			result = str(chrom) + "\t" + str(position)
			ref_base = position_data[2]
			nref_base = position_data[3]
			qual = float(position_data[4])
			genoKey = position_data[5]
			geno = position_data[6]
			#print base,
			if qual >= vcf_quality and qual != None:
				_p = sample.pileup(reference=chrom, start=position, end=position + 1)
				bases = list()
				for p in _p:
					if p.reference_pos == position:
						for r in p.pileups:
							if not r.is_del and not r.is_refskip:
								base = r.alignment.query_sequence[r.query_position-1]
								mapq = r.alignment.mapping_quality
								baseq = r.alignment.query_qualities[r.query_position-1]
								if mapq >= map_quality and baseq >= base_quality:
									bases.append(base)
			
				ref_count = 0
				depth = 0   
				for base in bases:
					depth += 1
					if base == ref_base:
						ref_count += 1
				alt_count = depth - ref_count
				gt_ind = genoKey.index("GT")
				ps_ind = None
				try:
					ps_ind = genoKey.index("PS")
				except ValueError:
					sys.stdout.write('%s is not in a phase set.' % result)
				if ps_ind is None:
					ps_val = "NA"
				else:
					ps_val = str(geno[ps_ind])
				#print str(depth) + "\t" + str(ref_count) + "\t" + str(alt_count)
				#if (int(ref_count) + int(alt_count)) > 0:
				result += "\t" + ref_base + "\t" + str(ref_count) + "\t" + nref_base + "\t" + str(alt_count) + '\t' + str(qual) + '\t' + str(geno[gt_ind]) + '\t' + ps_val
				print result
				i += 1
						
        
