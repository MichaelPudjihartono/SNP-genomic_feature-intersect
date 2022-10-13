import pandas as pd
import os
import sys
import pybedtools

import logging
import datetime
import argparse
import sys
import time

def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter SNPs based on selected features in BED format')        
    parser.add_argument(
        '-s', '--snps',
        help='''Input SNP file with the first three columns in BED format''' )
    parser.add_argument(
        '-sc', '--snps-codes3d',
        help='''Use this if you want to filter CoDeS3D output (i.e. significant_eqtls.txt) based on selected features''' )
    parser.add_argument(
        '-f', '--feature-dir', required=True,
        help='''Path to genomic feature directory, where bed files of each feature is stored in a named subdirectory, storing multiple files allowed''')
    parser.add_argument(
        '-o', '--output-dir', required=True,
        help='Directory to write results.')
    return parser.parse_args()
    


class Logger(object):
    def __init__(self, logfile=None, verbose=True):
        """
        A simple logger.
        """
        self.console = sys.stdout
        self.verbose = verbose
        if logfile is not None:
            self.log = open(logfile, 'w')
        else:
            self.log = None

        time = datetime.datetime.now()
        self.write((f'{time.strftime("%A")},'
                    f' {time.strftime("%b")}' 
                    f' {time.strftime("%d")},'
                    f' {time.strftime("%Y")}' 
                    f' {time.strftime("%I")}:'
                    f'{time.strftime("%M")}'
                    f' {time.strftime("%p")}'
        ))

    def write(self, message):
        if self.verbose:
            self.console.write(message+'\n')
        if self.log is not None:
            self.log.write(message+'\n')
            self.log.flush()

    def verbose(self, verbose=True):
        self.verbose = verbose    


#Called within the execution part as such:
#snp_bed = parse_codes3d(args.snps_codes3d, logger)
def parse_codes3d(snps_codes3d_args, logger):
    logger.write('Parsing codes3d output as BedTool object...')
    df = pd.DataFrame()
    
    if os.path.isfile(snps_codes3d_args):
        df = pd.read_csv(snps_codes3d_args, sep='\t')
        
        df = df[['snp_chr','snp_locus']].drop_duplicates()
        df = df.rename(columns={'snp_chr':'chrom', 'snp_locus':'chromEnd'})
        df['chromStart'] = df['chromEnd'].sub(1)
        df = df[['chrom','chromStart','chromEnd']]
        
        
        #convert df (aka codes3d output into BedTool object)
        snp_bed = pybedtools.BedTool.from_dataframe(df)

        return snp_bed

def parse_snps(snps_args, logger):
    logger.write('Parsing snps in BED format as BedTool object...')
    df = pd.DataFrame()
        


def parse_feature_df_codes3d(snps_codes3d_args):
    df = pd.read_csv(snps_codes3d_args, sep='\t') #There's no need to check if snps_codes3d_args is a file here since it's been checked above I think!
    
    feature_df = df[['snp','snp_chr', 'snp_locus']].drop_duplicates()
    feature_df = feature_df.rename(columns={'snp_chr':'chrom', 'snp_locus':'pos'})
    
    return feature_df
    
    
#Called below as such:
#feature(args.snps_codes3d, args.feature_dir, args.output_dir, logger)
def feature(snps_codes3d_args, feature_dir_args, output_dir_args, logger):

    snps_that_intersect_raw_features = parse_feature_df_codes3d(snps_codes3d_args)
    
    for feature in [d for d in os.listdir(feature_dir_args) if os.path.isdir(os.path.join(feature_dir_args, d))]: #This basically lists all file (aka "d") inside feature_dir_args that is a directory in itself
        x = feature.replace('COMPLEX_', '')
        logger.write(f'Identifying SNP overlap with {x}...')
        bedtool_instance_list = []
        
        if feature.startswith('COMPLEX_'):
            snps_that_intersect_all_subfeatures = parse_feature_df_codes3d(snps_codes3d_args)
            
            for subfeature in [d for d in os.listdir(os.path.join(feature_dir_args, feature)) if os.path.isdir(os.path.join(feature_dir_args, feature, d))]: #This lists all file (aka "d") inside the COMPLEX_ directory that iself is also a directory, e.g. [H3K4me3, H3K27ac]
                bedtool_subfeature_instance_list = []
                
                for subfeature_file in os.listdir(os.path.join(feature_dir_args, feature, subfeature)): #You already know that e.g. H3K27ac is a directory, now you iterate over the bedfiles in that directory
                    f = os.path.join(feature_dir_args, feature, subfeature, subfeature_file)
                    bedtool_subfeature_instance_list.append(pybedtools.BedTool(f))
                    
                snps_that_intersect_subfeature = snp_bed.intersect(bedtool_subfeature_instance_list, u=True).saveas(os.path.join(output_dir_args, f'snps_that_intersects_{subfeature}_subfeature_of_{feature}.bed'))
                
                if not os.path.getsize(os.path.join(output_dir_args, f'snps_that_intersects_{subfeature}_subfeature_of_{feature}.bed')) == 0:
                    snps_that_intersect_subfeature = pd.read_csv(os.path.join(output_dir_args, f'snps_that_intersects_{subfeature}_subfeature_of_{feature}.bed'), sep="\t", header = None)
                elif os.path.getsize(os.path.join(output_dir_args, f'snps_that_intersects_{subfeature}_subfeature_of_{feature}.bed')) == 0:
                    snps_that_intersect_subfeature = pd.DataFrame({0:[],1:[],2:[]})
                
                os.remove(os.path.join(output_dir_args, f'snps_that_intersects_{subfeature}_subfeature_of_{feature}.bed'))
                
                snps_that_intersect_subfeature[f'{subfeature}'] = True
                snps_that_intersect_subfeature.drop(1, inplace=True, axis=1) #drop the column corresponding to 'chromStart' in BED
                snps_that_intersect_subfeature.rename(columns={0:'chrom',2:'pos'}, inplace=True)
                
                snps_that_intersect_all_subfeatures = pd.merge(snps_that_intersect_all_subfeatures, snps_that_intersect_subfeature,  how='left', left_on=['chrom','pos'], right_on = ['chrom','pos']).fillna(False)
                
            snps_that_intersect_all_subfeatures = snps_that_intersect_all_subfeatures[(snps_that_intersect_all_subfeatures.iloc[:, 3:] == True).all(axis=1)]
            
            
            snps_that_intersect_all_subfeatures[f'{x}'] = True
            snps_that_intersect_all_subfeatures = snps_that_intersect_all_subfeatures[['chrom', 'pos', f'{x}']]
            
            snps_that_intersect_raw_features = pd.merge(snps_that_intersect_raw_features, snps_that_intersect_all_subfeatures,  how='left', left_on=['chrom','pos'], right_on = ['chrom','pos']).fillna(False)
            
            
        elif not feature.startswith('COMPLEX_'):
            for feature_file in os.listdir(os.path.join(feature_dir_args, feature)):
                f = os.path.join(feature_dir_args, feature, feature_file)
                bedtool_instance_list.append(pybedtools.BedTool(f))
                
            snps_that_intersect_feature = snp_bed.intersect(bedtool_instance_list, u=True).saveas(os.path.join(output_dir_args, f'snps_that_intersect_{feature}.bed'))
            
            if not os.path.getsize(os.path.join(output_dir_args, f'snps_that_intersect_{feature}.bed')) == 0:
                snps_that_intersect_feature = pd.read_csv(os.path.join(output_dir_args, f'snps_that_intersect_{feature}.bed'), sep="\t", header = None)
            elif os.path.getsize(os.path.join(output_dir_args, f'snps_that_intersect_{feature}.bed')) == 0:
                snps_that_intersect_feature = pd.DataFrame({0:[],1:[],2:[]})
            
            os.remove(os.path.join(output_dir_args, f'snps_that_intersect_{feature}.bed'))
            
            snps_that_intersect_feature[f'{feature}'] = True
            snps_that_intersect_feature.drop(1, inplace=True, axis=1) #drop the column corresponding to 'chromStart' in BED
            snps_that_intersect_feature.rename(columns={0:'chrom',2:'pos'}, inplace=True)
            
            snps_that_intersect_raw_features = pd.merge(snps_that_intersect_raw_features, snps_that_intersect_feature,  how='left', left_on=['chrom','pos'], right_on = ['chrom','pos']).fillna(False)
            
    return snps_that_intersect_raw_features
            
    


if __name__=='__main__':
    args = parse_args()
    if not args.snps and not args.snps_codes3d: #This checks if you have at least one of the required input types
        sys.exit('FATAL: One of --snps or --snps-codes3d is required.\nExiting.')
  
    start_time = time.time()
    
    if args.snps and args.snps_codes3d: #Another check! have to have at least one of the input type
        sys.exit('Only one of --snps or --snps-codes3d is required.\nExiting.')

    os.makedirs(args.output_dir, exist_ok=True) #This makes a new directory as specified in the -o argument #exist_ok parameter make sure that if the output_dir already exist, no error will be raised #If the target directory already exists an OSError is raised if its value is False otherwise not. For value True leaves directory unaltered.
    
    logger = Logger(logfile=os.path.join(args.output_dir, 'feature.log')) #This instantiate the Logger class as defined above, with parameter logfile specified

    logger.write('SETTINGS\n========')
    for arg in vars(args): #args = parse_args(), 
    #This iterates over all arguments called when you execute the code in linux
    #The vars() function returns the __dict__ attribute of the given object. So it iterates over the args object (which was previously defined as args = parse_args())
        logger.write(f'{arg}:\t {getattr(args, arg)}')
    logger.write('\n')
    
    if args.snps_codes3d:
        snp_bed = parse_codes3d(args.snps_codes3d, logger)
        
        snps_that_intersect_raw_features = feature(args.snps_codes3d, args.feature_dir, args.output_dir, logger)
        snps_that_intersect_raw_features.to_csv(os.path.join(args.output_dir, 'raw_feature_intersect.txt'), sep='\t', index=False)
        
        #Output only SNPs that pass filter as a file, one rsID each row
        snps_that_pass_filter = snps_that_intersect_raw_features[(snps_that_intersect_raw_features.iloc[:, 3:] == True).any(axis=1)]
        snps_that_pass_filter = snps_that_pass_filter['snp']
        snps_that_pass_filter.to_csv(os.path.join(args.output_dir, 'snps_that_pass_filter.txt'), index=False, header=False)
        
        #Create a filtered version of significant_eqtls.txt
        logger.write(f'Writing FILTERED codes3d output...')
        snps_that_pass_filter_set = set(snps_that_pass_filter)
        df = pd.read_csv(args.snps_codes3d, sep="\t")
        df_filtered = df[df['snp'].isin(snps_that_pass_filter_set)]
        
        df_filtered.to_csv(os.path.join(args.output_dir, 'significant_eqtls.txt'), sep='\t', index=False)
    
    if args.snps:
        snp_bed = parse_snps(args.snps, logger)
        
    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60: .2f} minutes.')
        
        
        
    
    
    
    
    
    
    
    







    
    
    
    
    
    
    
    
 #multimorbid3D:
    #The concatenated dataframe (aka "df") outputed from the LD expansion have no row of the query_snps itself as the target
    #This is an algorithm to make a (currently separate) dataframe of each query_snp that exist in the ld db as both the query and target, with each of the ['chromq',	'posq', 'rsidq', 'chromt', 'post', 'rsidt', 'corr','dprime'] properly filled
    #query_snps_df = pd.DataFrame({'rsidq': query_snps,
                        #'corr': 1,
                        #'dprime': 1
                       #})                   
    #query_snps_df = pd.merge(df[['chromq','posq','rsidq']], query_snps_df, on='rsidq').drop_duplicates()
    #query_snps_df[['chromt','rsidt','post']]=query_snps_df[['chromq','rsidq','posq']]
    #query_snps_df = query_snps_df[['chromq', 'posq', 'rsidq', 'chromt', 'post', 'rsidt', 'corr', 'dprime']]
    
                     
    #This concatenates the dataframe of the original query snps and "df" (aka the concatenated dataframe output of the LD expansion) and sort the total giant df                   
    #df = (pd.concat([df, query_snps_df]) 
             #.sort_values(by=['rsidq', 'corr', 'dprime'], ascending=False))