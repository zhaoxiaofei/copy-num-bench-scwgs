import argparse
import logging
import sys

logging.basicConfig(level=logging.DEBUG)

def myint(n): return int(round(float(n)))

def main():
    prefixCN = ''
    parser = argparse.ArgumentParser(description='Transform raw single-cell copy-number calling result into the standard BED format. ')
    parser.add_argument('--caller', type=str, required=True, help='Single-cell copy-number caller (can be set to hmmcopy, ginkgo, etc.). ')
    parser.add_argument('--sample', type=str, required=True, help='Keyword that can uniquely identify a sample (for example, SRR927014). ')
    parser.add_argument('--dp', action='store_true', help='Keyword that can uniquely identify a sample (for example, SRR927014). ')
    parser.add_argument('--fai', type=str, default='', help='The fai (fasta index) file of the reference')

    # parser.add_argument('--infile', type=str, required=True, help='The input file. ') 
    args = parser.parse_args()
    
    chrom_sizes = {}
    if args.fai:
        with open(args.fai) as file:
            for line in file:
                toks = line.strip().split()
                chrom_name = str(toks[0].strip('"').strip("'"))
                chrom_size = int(toks[1].strip('"').strip("'"))
                chrom_sizes[chrom_name] = chrom_size
    bed_header = '\t'.join(['#chr_37', 'start_37', 'end_37', ('obsDP' if args.dp else 'obsCN')])
    print(bed_header)
    
    if   args.caller == 'hmmcopy': # 0/segs.csv for both int CN and non-int CN
        for line in sys.stdin:
            tokens = line.strip().split(',')
            if tokens[0] == 'chr':
                continue
            chrom = tokens[0]
            start = myint(tokens[1])
            end = myint(tokens[2])
            cn = float(tokens[3]) if args.dp else myint(tokens[3])
            print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn}')
    elif args.caller == 'ginkgo': # SegCopy (int CN) or SegFixed (non-int CN)
        sample_col_index = 0
        for line in sys.stdin:
            tokens = line.strip().split()
            if tokens[0] == 'CHR':
                for index, token in enumerate(tokens):
                    if args.sample in token: sample_col_index = index
                assert sample_col_index > 2, F'The sample keyword {args.sample} is not found!'
                continue
            chrom = tokens[0]
            start = myint(tokens[1])
            end = myint(tokens[2])
            cn = float(tokens[sample_col_index]) if args.dp else myint(tokens[sample_col_index])
            print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn}')
    elif args.caller == 'copynumber': # output.csv
        # CopyNumber does not provide integer copy numbers, it only provides relative signal intensity. 
        # Therefore, we assume that the last-column sample is diploid normal and estimate copy number this way. 
        sample_col_index = 0
        for line in sys.stdin:
            tokens = line.strip().split(',')
            if tokens[0] == 'chrom' or tokens[0] == '"chrom"':
                for index, token in enumerate(tokens):
                    if args.sample in token: sample_col_index = index
                assert sample_col_index > 2, F'The sample keyword {args.sample} is not found!'
                continue
            assert sample_col_index > 0, F'The sample keyword {args.sample} is not found!'
            assert len(tokens) > 5, 'At least two cells are required for copynumber!'
            chrom_index = int(tokens[0])
            chrom = F'chr{chrom_index}'
            start = myint(tokens[2])
            end = myint(tokens[3])
            if args.dp:
                cn_float = float(tokens[sample_col_index])
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{(cn_float)}')
            else:
                normal_diploid_index = ((len(tokens)-2) if (len(tokens) - 1 == sample_col_index) else (len(tokens)-1))
                cn_float = (float(tokens[sample_col_index]) + 1.0) / (float(tokens[normal_diploid_index]) + 1.0) * 2
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{round(cn_float)}')
    elif args.caller == 'sccnv': 
        ''' 
        File: result.dat6_cnvsmooth.txt
        Content:
            chr	pos1	pos2	Mappability	AT	GC	N	\
                SRR927014_12.sort-rmdup.pb2.simulating-HCC1395.mapq30	SRR927015_12.sort-rmdup.fpn.simulating-HCC1395.mapq30	\
                SRR927017_12.sort-rmdup.pb2.simulating-HeLa.mapq30	SRR927018_12.sort-rmdup.fpn.simulating-HeLa.mapq30	\
                SRR927019_12.sort-rmdup.cn2.simulating-normal.mapq30	SRR927025_12.sort-rmdup.pb2.simulating-HCC1395.mapq30	SRR927026_12.sort-rmdup.fpn.simulating-HCC1395.mapq30
            chr1	500000	1000000	0.545064	0.47728	0.479986	0.042736	\
                1.9930497674358882	1.902179715554554	2.1670958230016883	2.220547643581603	2.0684420508394425	1.9005440903049953	2.02808425386971
            etc. 
        '''
        sample_col_index = 0
        for line in sys.stdin:
            tokens = line.strip().split()
            if tokens[0] == 'chr':
                for index, token in enumerate(tokens):
                    if args.sample in token: sample_col_index = index
                assert sample_col_index > 2, F'The sample keyword {args.sample} is not found!'
                continue
            chrom = tokens[0]
            start = myint(tokens[1])
            end = myint(tokens[2])
            cn_string = tokens[sample_col_index]
            if 'nan' != cn_string:
                if args.dp:
                    print(F'{chrom}\t{start}\t{end}\t{prefixCN}{(float(cn_string))}')
                else:
                    print(F'{chrom}\t{start}\t{end}\t{prefixCN}{round(float(cn_string))}')
    elif args.caller == 'scnv':
        # SCNV has no doc, so we did not run it yet. 
        pass
    elif args.caller == 'secnv':
        if args.dp:
            '''
            File: SRA091188.out.dir/S01.secnv.dir/output/genome_cov.bed
            Content:
            chromosome      start   stop    SRR926893_12    SRR926893_12
            chr1    1000001 1500000 392     707     536     809     261     450
            '''
            sample_col_index = 0
            for line in sys.stdin:
                tokens = line.strip().split()
                if tokens[0] == 'chromosome':
                    for index, token in enumerate(tokens):
                        if args.sample in token: sample_col_index = index
                    assert sample_col_index > 2, F'The sample keyword {args.sample} is not found!'
                    continue
                chrom = tokens[0]
                start = myint(tokens[1])
                end = myint(tokens[2])
                cn = float(tokens[sample_col_index]) # if args.dp else myint(tokens[sample_col_index])
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn}')
        else:
            '''
            File: cnv_matrix.csv
            Content: 
              ,chr1:1000001-1500000,chr1:1500001-2000000,etc.
              SRR927014_12,4.0,4.0,etc.
              etc. 
            '''
            row2sample_col2region = [line.strip().split(',') for line in sys.stdin]
            sample_row_index = -1
            for i, sample2regions in enumerate(row2sample_col2region):
                sample = sample2regions[0]
                if args.sample in sample: sample_row_index = i 
            assert sample_row_index > 0, F'The sample keyword {args.sample} is not found!'
            for region, cn in zip(row2sample_col2region[0], row2sample_col2region[sample_row_index]):
                if region.strip() == '': continue
                chrom = region.split(':')[0]
                start = myint(region.split(':')[1].split('-')[0])
                end = myint(region.split(':')[1].split('-')[1])
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{round(float(cn))}')
    elif args.caller == 'scope' or args.caller == 'scyn':
        '''
        File: scyn_output.csv
        Content: 
          ,SRR927014_12,SRR927015_12,SRR927017_12,SRR927018_12,SRR927019_12,SRR927024_12
          chr1:1000001-1500000,2,2,3,3,2,2
          etc.
        '''
        sample_col_index = 0
        for line in sys.stdin:
            tokens = line.strip().split(',')
            if tokens[0] == '':
                for index, token in enumerate(tokens):
                    if args.sample in token: sample_col_index = index
                assert sample_col_index > 0, F'The sample keyword {args.sample} is not found!'
                continue
            chrom = tokens[0].split(':')[0]
            start = myint(tokens[0].split(':')[1].split('-')[0])
            end = myint(tokens[0].split(':')[1].split('-')[1])
            cn_string = tokens[sample_col_index]
            if args.dp: # SCYN always generate integer CN
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn_string}')
            else:
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn_string}')
    elif args.caller == 'sconce2':
        for line in sys.stdin:
            tokens = line.strip().split()
            if args.dp:
                print(F'{tokens[0]}\t{tokens[1]}\t{tokens[2]}\t{prefixCN}{(float(tokens[3]))}')
            else:
                print(F'{tokens[0]}\t{tokens[1]}\t{tokens[2]}\t{prefixCN}{round(float(tokens[3]))}')
    elif args.caller == 'chisel':
        mode = None
        sample_barcode = None
        for line in sys.stdin:
            if line.startswith('#CELL'): 
                mode = 'barcode'
                continue
            if line.startswith('#CHR'): 
                mode = 'cnv'
                continue
            tokens = line.strip().split()
            if mode == 'barcode':
                if args.sample in tokens[0]: sample_barcode = tokens[1]
            if mode == 'cnv' and tokens[3] == sample_barcode:
                RDR = float(tokens[6])
                majorCN, minorCN = tuple(int(n) for n in tokens[12].split('|'))
                if args.dp:
                    print(F'{tokens[0]}\t{tokens[1]}\t{tokens[2]}\t{prefixCN}{RDR}')
                else:
                    print(F'{tokens[0]}\t{tokens[1]}\t{tokens[2]}\t{prefixCN}{majorCN+minorCN}')
    elif args.caller == 'alleloscope':
        '''
        File: allelo_out.csv  (output of simplerun_alleloscope_v21.R)
        Content (CSV with header):
          cell,chrN:start-end_rho_hat,...,chrN:start-end_theta_hat,...,chrN:start-end_genotype_state,...
          cellname1,1.02,...,0.51,...,3,...

        rho_hat ≈ total_CN / 2  (1.0 = diploid)
        genotype_state: 1=homdel(CN0), 2=LOH(CN1), 3=diploid(CN2),
                        4=cnLOH(CN2), 5=gain(CN3), 6=amp(CN4)
        For --dp mode: output rho_hat * 2 (continuous CN estimate)
        For int mode:  output round(rho_hat * 2)
        '''
        import csv as csv_mod
        reader = csv_mod.DictReader(sys.stdin)
        # Find the row matching --sample
        for row in reader:
            cell = row.get('cell', '')
            if args.sample not in cell:
                continue
            # Extract rho_hat columns (they contain the CN signal)
            for col_name, col_val in row.items():
                if '_rho_hat' not in col_name:
                    continue
                # Column name format: "chrN:start-end_rho_hat"
                region_str = col_name.replace('_rho_hat', '')
                chrom = region_str.split(':')[0]
                coords = region_str.split(':')[1]
                start = myint(coords.split('-')[0])
                end   = myint(coords.split('-')[1])
                rho = float(col_val)
                cn_continuous = rho * 2.0
                if args.dp:
                    print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn_continuous}')
                else:
                    print(F'{chrom}\t{start}\t{end}\t{prefixCN}{round(cn_continuous)}')
            break  # only one matching row needed
    elif args.caller == 'aneufinder':
        # Decompressed BED from AneuFinder: multiple tracks, each with lines like:
        # chr1	0	26000000	3-somy	0	.	0	26000000	205,0,0
        found_sample = False
        for line in sys.stdin:
            if line.startswith('track '):
                found_sample = (args.sample in line)
                continue
            if not found_sample:
                continue
            tokens = line.strip().split('\t')
            if len(tokens) < 4:
                continue
            chrom = tokens[0]
            start = myint(tokens[1])
            end   = myint(tokens[2])
            somy_str = tokens[3]  # e.g. "3-somy"
            cn = float(somy_str.split('-')[0])
            if args.dp:
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{cn}')
            else:
                print(F'{chrom}\t{start}\t{end}\t{prefixCN}{myint(cn)}')
    elif args.caller == 'flcna':
        # FLCNA output CSV: sampleID,samplename,Cluster,chr,start,end,state,start.coor,end.coor,width_bins
        # States: Del.d=0, Del.s=1, normal=2 (not listed), Dup.s=3, Dup.d=4
        state2cn = {'Del.d': 0, 'Del.s': 1, 'Dup.s': 3, 'Dup.d': 4}
        segments = []  # collect all aberrant segments first
        for line in sys.stdin:
            tokens = line.strip().split(',') if ',' in line else line.strip().split()
            tokens = [t.strip('"') for t in tokens]
            if len(tokens) < 9:
                continue
            if tokens[0] == 'sampleID' or tokens[1] == 'samplename':
                continue
            samplename = tokens[1]
            if args.sample not in samplename:
                continue
            chrom = tokens[3].strip().strip('"').strip("'")
            start = myint(tokens[7])
            end   = myint(tokens[8])
            state = tokens[6]
            cn = state2cn.get(state, 2)
            if start < end:
                segments.append((chrom, start, end, cn))
            else:
                logging.warning(f'The segment {(chrom, start, end, cn)} from sample={samplename} is invalid! Skipping this segment!')
        # Group by chromosome and sort
        from collections import defaultdict
        chrom_segs = defaultdict(list)
        for chrom, start, end, cn in segments:
            chrom_segs[chrom].append((start, end, cn))
        for chrom in chrom_segs:
            chrom_segs[chrom].sort()
            segs = [seg for seg in chrom_segs[chrom]]
            new_segs = []
            # ensure that BED regions are non-overlapping by trimming each region if needed
            for i in range(len(segs)):
                max_start = segs[i][0]
                min_end = segs[i][1]
                for j in range(0, i, 1): max_start = max((max_start, segs[j][1]))
                for j in range(i+1, len(segs), 1): min_end = min((min_end, segs[j][0]))
                if max_start < min_end:
                    new_segs.append((max_start, min_end, segs[i][2]))
            chrom_segs[chrom] = new_segs
        # https://sorryios.ai/chat/e0377335-c9a1-43ba-bb36-5143c3479f4a
        # Including normal copy number regions in FLCNA benchmarking
        # Fill gaps with CN=2
        for chrom, segs in chrom_segs.items(): # sorted(chrom_segs.items(), key=lambda x: x[0]):
            prev_end = 0  # or use 1 if 1-based
            for start, end, cn in segs:
                if start > prev_end:
                    # gap → normal CN=2
                    if args.dp:
                        print(f'{chrom}\t{prev_end}\t{start}\t{prefixCN}{2.0}')
                    else:
                        print(f'{chrom}\t{prev_end}\t{start}\t{prefixCN}2')
                if args.dp:
                    print(f'{chrom}\t{start}\t{end}\t{prefixCN}{float(cn)}')
                else:
                    print(f'{chrom}\t{start}\t{end}\t{prefixCN}{cn}')
                prev_end = end
            # Optionally: tail region after last segment to chrom end
            # if you have chrom_sizes dict:
            if chrom_sizes:
                assert chrom in chrom_sizes, f'The chromosome {chrom} is not in {chrom_sizes}'
            if chrom_sizes and prev_end < chrom_sizes[chrom]:
                print(f'{chrom}\t{prev_end}\t{chrom_sizes[chrom]}\t{prefixCN}2')
    else:
        logging.fatal(F'The copy-number caller {args.caller} is invalid!')

if __name__ == '__main__': main()
