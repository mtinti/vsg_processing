#!/usr/bin/env python3
"""
VSG Meta-Analysis Pipeline
===========================
A comprehensive pipeline for analyzing VSG expression data from RNA-seq experiments.

This script processes RNA-seq count data to:
1. Generate experiment configuration (exp_config.json)
2. Calculate silent VSG expression changes (silentCsvData.csv)
3. Calculate main VSG expression changes (mainCsvData.csv)
4. Generate quality control metrics (QC.csv)

Author: [Your Name]
Date: 2025
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Suppress pandas warnings
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)

def make_desc(_GFF):
    gff =pd.read_csv( _GFF, sep='\t', header=None, comment='#')
    #print(gff)
    gff = gff[gff.iloc[:,2].str.contains('gene')]
    #print( gff[gff[gff.columns[-1]].str.contains('Tb427_020006200')] )
    #print(gff)
    desc = {}
    chr_dict = {}
    start_dict = {}
    sense_dict = {}
    for index,n in enumerate(gff.iloc[:,-1]):
        n=n.replace('%2C',' ')
        item_list = n.split(';')
        #print (item_list)
        temp_dict = {}
        for m in item_list:
            #print(m)
            temp_dict[m.split('=')[0].strip()]=m.split('=')[1].strip()
        #print(temp_dict['ID'])
        #print(temp_dict['description'])
        desc[temp_dict['ID']]=temp_dict.get('description','none')
        chr_dict[temp_dict['ID']]=gff.iloc[index][0]
        start_dict[temp_dict['ID']]=gff.iloc[index][3]
        sense_dict[temp_dict['ID']]=gff.iloc[index][6]
        
    return desc,chr_dict,start_dict,sense_dict

desc_dict_gff, _, _,_ = make_desc('TriTrypDB-68_TbruceiTREU927.gff')
#print(desc_dict_gff['Tb927.1.4900'])



class VSGAnalysisPipeline:
    """Main pipeline class for VSG expression analysis."""
    
    def __init__(self, 
                 experiment_table_path: str = 'experiment_table.txt',
                 vsg_dict_path: str = 'vsg_dic.txt',
                 base_result_path: str = 'myRna-seq/results/result_vsgs',
                 output_dir: str = 'output'):
        """
        Initialize the VSG analysis pipeline.
        
        Args:
            experiment_table_path: Path to experiment metadata table
            vsg_dict_path: Path to VSG dictionary file
            base_result_path: Base path for result files
            output_dir: Directory for output files
        """
        self.experiment_table_path = experiment_table_path
        self.vsg_dict_path = vsg_dict_path
        self.base_result_path = base_result_path
        self.output_dir = Path(output_dir)
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.df = None
        self.vsg_dict = None
        self.exp_config = None
        
        # Samples to exclude (known issues)
        self.exclude_samples = {
            'SRR7401110', 'SRR7401109',  # PIP5Pase_1: mdsum could not be verified
            'ERR4160685', 'ERR4160687', 'ERR4160684', 'ERR4160686'  # ZC3H5_1: control not deposited
        }
    
    def load_experiment_table(self) -> pd.DataFrame:
        """Load and preprocess the experiment table."""
        logger.info("Loading experiment table...")
        
        df = pd.read_csv(self.experiment_table_path, sep='\t')
        df = df.dropna(subset=['FASTQ_TYPE'])
        
        # Map FASTQ types
        tmp = df['FASTQ_TYPE'].copy()
        df['FASTQ_TYPE'] = df['FASTQ_TYPE'].map({
            'PAIRED': 'paired',
            'SINGLE': 'single',
            'NANOPORE': 'single'
        })
        df['READ_TYPE'] = tmp.map({
            'PAIRED': 'paired',
            'SINGLE': 'single',
            'NANOPORE': 'nanopore'
        })
        
        # Add file paths
        df['count_file'] = df.apply(
            lambda row: f"{self.base_result_path}/{row['SRA_id']}/{row['SRA_id']}_counts_{row['FASTQ_TYPE']}_all.txt",
            axis=1
        )
        df['qc_file'] = df.apply(
            lambda row: f"{self.base_result_path}/{row['SRA_id']}/qc/feature_counts/{row['SRA_id']}_counts_{row['FASTQ_TYPE']}_all.txt.summary",
            axis=1
        )
        
        self.df = df
        logger.info(f"Loaded {len(df)} samples")
        return df
    
    def load_vsg_dictionary(self) -> Dict[str, str]:
        """Load VSG dictionary mapping chromosomes to families."""
        logger.info("Loading VSG dictionary...")
        
        vsg_dict = pd.read_csv(self.vsg_dict_path, sep=' ', header=None)
        vsg_dict.columns = ['Chr', 'Family']
        vsg_dict['Chr'] = vsg_dict['Chr'].str.replace('>', '')
        vsg_dict['Family'] = vsg_dict['Family'].str.split('_').str[0]
        
        self.vsg_dict = vsg_dict.set_index('Chr').to_dict()['Family']
        logger.info(f"Loaded {len(self.vsg_dict)} VSG mappings")
        return self.vsg_dict
    
    def validate_files(self) -> None:
        """Validate that all required files exist."""
        logger.info("Validating files...")
        
        missing_count_files = []
        missing_qc_files = []
        
        for idx, row in self.df.iterrows():
            if not Path(row['count_file']).exists():
                missing_count_files.append(row['count_file'])
            if not Path(row['qc_file']).exists():
                missing_qc_files.append(row['qc_file'])
        
        if missing_count_files:
            logger.warning(f"Missing {len(missing_count_files)} count files")
        if missing_qc_files:
            logger.warning(f"Missing {len(missing_qc_files)} QC files")
    
    def get_count_file(self, file_path: str, filter_vsg: bool = False) -> pd.DataFrame:
        """
        Load and process a count file.
        
        Args:
            file_path: Path to count file
            filter_vsg: Whether to filter for VSG genes only
        
        Returns:
            DataFrame with count data
        """
        tmp = pd.read_csv(file_path, sep='\t', comment='#')
        
        if filter_vsg:
            tmp = tmp[tmp['Chr'].str.startswith('Tb427VSG-')]
        
        cols = list(tmp.columns)
        cols[-1] = 'counts'
        tmp.columns = cols
        tmp['Family'] = tmp['Chr'].map(self.vsg_dict)
        
        return tmp
    
    def get_assigned_reads(self, file_path: str) -> int:
        """Get number of assigned reads from QC file."""
        tmp = pd.read_csv(file_path, sep='\t')
        return tmp[tmp['Status'] == 'Assigned'].iloc[:, -1].values[0]
    
    def make_relative_counts(self, count_path: str, qc_path: str) -> pd.Series:
        """Calculate relative counts (counts/total assigned reads)."""
        counts = self.get_count_file(count_path)
        total =  counts[ counts['Chr'].isin( [
            'Tb927_01_v5.1',
            'Tb927_02_v5.1',
            'Tb927_03_v5.1',
            'Tb927_04_v5.1',
            'Tb927_05_v5.1',
            'Tb927_06_v5.1',
            'Tb927_07_v5.1',
            'Tb927_08_v5.1',
            'Tb927_09_v5.1',
            'Tb927_10_v5.1',
            'Tb927_11_v5.1',
            
        ]
                                           ) ]['counts'].sum()  #self.get_assigned_reads(qc_path)
        return counts['counts'] / total
    
    def identify_top_vsgs(self) -> None:
        """Identify the top expressed VSG for each sample."""
        logger.info("Identifying top VSGs...")
        
        top_vsgs = []
        qc_vsgs = []
        
        for _, row in tqdm(self.df.iterrows(), total=len(self.df), desc="Processing samples"):
            counts = self.get_count_file(row['count_file'], filter_vsg=True)
            
            # Get top VSG
            top_vsg = counts.nlargest(1, 'counts')['Chr'].values[0]
            top_vsgs.append(top_vsg)
            
            # Calculate QC metric (ratio of second to top VSG)
            top_two = counts.nlargest(2, 'counts')['counts'].values
            if len(top_two) >= 2:
                qc_ratio = top_two[1] / (top_two[0] + top_two[1])
            else:
                qc_ratio = 0
            qc_vsgs.append(qc_ratio)
        
        self.df['top_vsgs'] = top_vsgs
        self.df['qc_vsgs'] = qc_vsgs
        
        logger.info(f"Top VSG distribution:\n{self.df['top_vsgs'].value_counts().head()}")
    
    def generate_experiment_config(self) -> Dict:
        """Generate experiment configuration dictionary."""
        logger.info("Generating experiment configuration...")
        
        # Filter out PCF samples and excluded samples
        df_filtered = self.df[
            (self.df['notes'] != 'PCF') & 
            (~self.df['SRA_id'].isin(self.exclude_samples))
        ]
        
        exp_config = {}
        
        for exp_name, exp_df in df_filtered.groupby('Proposed_Experiment_Name'):
            # Find control samples
            control_df = exp_df[exp_df['Proposedcondition'].isin(['Control', 'ControlUN'])]
            
            if control_df.empty:
                logger.warning(f"No controls found for {exp_name}")
                continue
            
            # Get metadata from controls
            top_vsg = control_df['top_vsgs'].values[0]
            pubmed_id = control_df['PMID'].values[0]
            title = control_df['Title'].values[0]
            fastq_type = control_df['READ_TYPE'].values[0]
            
            # Process each treatment condition
            treatment_df = exp_df[~exp_df['Proposedcondition'].isin(['Control', 'ControlUN'])]
            
            for condition, cond_df in treatment_df.groupby('Proposedcondition'):
                exp_key = f"{exp_name}_{condition}"
                exp_config[exp_key] = {
                    'experiment':exp_key,
                    'controls': control_df['SRA_id'].tolist(),
                    'treatments': cond_df['SRA_id'].tolist(),
                    'top_vsg': [top_vsg],
                    'pubmed_id': int(pubmed_id),
                    'title': title,
                    'fastq_type': fastq_type
                }
        
        # Apply manual corrections for known mixed populations
        if 'PIP5Pase_1_KO' in exp_config:
            exp_config['PIP5Pase_1_KO']['top_vsg'] = ['Tb427VSG-2', 'Tb427VSG-11']
        if 'PIP5Pase_2_D360A_N362A' in exp_config:
            exp_config['PIP5Pase_2_D360A_N362A']['top_vsg'] = ['Tb427VSG-2', 'Tb427VSG-567']
        
        self.exp_config = exp_config
        logger.info(f"Generated configuration for {len(exp_config)} experiments")
        return exp_config
    
    def process_experiment(self, exp_dict: Dict) -> pd.DataFrame:
        """
        Process a single experiment to calculate log2 fold changes.
        
        Args:
            exp_dict: Experiment configuration dictionary
        
        Returns:
            DataFrame with log2 fold changes for each gene
        """
        # Initialize merged dataframes
        template = self.get_count_file(self.df['count_file'].iloc[0])
        merged_control = template.copy()
        merged_control['counts'] = 0
        merged_treatment = merged_control.copy()
        
        # Sum control samples
        for sra_id in exp_dict['controls']:
            row = self.df[self.df['SRA_id'] == sra_id].iloc[0]
            relative_counts = self.make_relative_counts(row['count_file'], row['qc_file'])
            merged_control['counts'] += relative_counts
        
        # Sum treatment samples
        for sra_id in exp_dict['treatments']:
            row = self.df[self.df['SRA_id'] == sra_id].iloc[0]
            relative_counts = self.make_relative_counts(row['count_file'], row['qc_file'])
            merged_treatment['counts'] += relative_counts
        
        # Average by number of replicates
        merged_control['counts'] /= len(exp_dict['controls'])
        merged_treatment['counts'] /= len(exp_dict['treatments'])
        
        # Calculate log2 fold change
        result = template.copy()
        #result['counts'] = (
        #    np.log2(merged_treatment['counts'].values * 1000000 + 1) - 
        #    np.log2(merged_control['counts'].values * 1000000 + 1)
        #)
        # Calculate fold change
        result['counts'] = (
            (merged_treatment['counts'].values * 1000000 + 1) /
            (merged_control['counts'].values * 1000000 + 1)
        ) 

        de_analysis = result.copy()
        de_analysis['log_2FC'] = (
            np.log2(merged_treatment['counts'].values * 1000000 + 1) - np.log2(merged_control['counts'].values * 1000000 + 1)
        ) 
        de_analysis['desc']=[desc_dict_gff.get(n,n) for n in de_analysis['Geneid']]
        
        fname = exp_dict['experiment']
        #print(exp_dict)
        de_analysis.to_csv(f'de_analysis/{fname}.csv')
        # Filter for VSGs only
        result = result[result['Chr'].str.startswith('Tb427VSG-')]
        
        return result
    
    def analyze_silent_vsgs(self) -> pd.DataFrame:
        """Analyze silent VSG expression changes across all experiments."""
        logger.info("Analyzing silent VSG expression...")
        
        collect_results = []
        
        for exp_name in tqdm(self.exp_config, desc="Processing experiments"):
            result = self.process_experiment(self.exp_config[exp_name])
            
            # Exclude main VSGs
            result = result[~result['Chr'].isin(self.exp_config[exp_name]['top_vsg'])]
            
            # Keep only upregulated VSGs
            result = result[result['counts'] > 0]
            
            # Sum by family
            #family_sums = result.groupby('Family')['counts'].sum()
            
            #family_sums = result.groupby('Family')['counts'].quantile(q=0.25)
            #family_sums = result.groupby('Family')['counts'].mean()
            family_sums = result.groupby('Family')['counts'].median()
            
            family_sums.name = exp_name
            collect_results.append(family_sums)
        
        # Combine results
        results_df = pd.concat(collect_results, axis=1).fillna(0).T
        results_df['sum'] = results_df.sum(axis=1)
        
        logger.info(f"Silent VSG analysis complete. Shape: {results_df.shape}")
        return results_df
    
    def analyze_main_vsgs(self) -> pd.DataFrame:
        """Analyze main VSG expression changes across all experiments."""
        logger.info("Analyzing main VSG expression...")
        
        collect_results = []
        
        for exp_name in tqdm(self.exp_config, desc="Processing experiments"):
            result = self.process_experiment(self.exp_config[exp_name])
            
            # Keep only main VSGs
            result = result[result['Chr'].isin([self.exp_config[exp_name]['top_vsg'][0]])]
            #get only first in case of two main vsgs (VSG2 x experiment PIP5Pase_2_D360A_N362A and PIP5Pase_1_KO)

            #if result.shape[0] >1:
            #    result=result.iloc[:0,:]
            #if exp_name == 'PIP5Pase_2_D360A_N362A':
            #    print(result)            
            
            total_change = np.log2(result[['counts']].mean())
            total_change.name = exp_name
            collect_results.append(total_change)
        
        # Combine results
        results_df = pd.concat(collect_results, axis=1).fillna(0).T
        
        logger.info(f"Main VSG analysis complete. Shape: {results_df.shape}")
        return results_df
    
    def calculate_qc_metrics(self) -> pd.DataFrame:
        """Calculate quality control metrics for each experiment."""
        logger.info("Calculating QC metrics...")
        
        qc_results = []
        
        for exp_name in tqdm(self.exp_config, desc="Calculating QC"):
            main_vsg_perc = self._calculate_main_vsg_percentage(self.exp_config[exp_name])
            mito_perc = self._calculate_mito_percentage(self.exp_config[exp_name])
            
            qc_results.append({
                'experiment': exp_name,
                'Main_VSG_perc': main_vsg_perc,
                'Mito_perc': mito_perc
            })
        
        qc_df = pd.DataFrame(qc_results)
        logger.info(f"QC metrics calculated for {len(qc_df)} experiments")
        return qc_df
    
    def _calculate_main_vsg_percentage(self, exp_dict: Dict) -> float:
        """Calculate percentage of reads mapping to main VSG."""
        percentages = []
        
        for sra_id in exp_dict['controls']:
            row = self.df[self.df['SRA_id'] == sra_id].iloc[0]
            counts = self.get_count_file(row['count_file'])
            relative_counts = self.make_relative_counts(row['count_file'], row['qc_file'])
            
            counts['relative_counts'] = relative_counts
            counts_indexed = counts.set_index('Geneid')
            
            main_vsg_genes = [f'gene_{vsg}' for vsg in exp_dict['top_vsg']]
            if all(gene in counts_indexed.index for gene in main_vsg_genes):
                total = counts_indexed.loc[main_vsg_genes]['relative_counts'].sum()
                percentages.append(total)
        
        return np.mean(percentages) if percentages else 0
    
    def _calculate_mito_percentage(self, exp_dict: Dict) -> float:
        """Calculate percentage of reads mapping to mitochondrial genes."""
        percentages = []
        
        for sra_id in exp_dict['controls']:
            row = self.df[self.df['SRA_id'] == sra_id].iloc[0]
            counts = self.get_count_file(row['count_file'])
            relative_counts = self.make_relative_counts(row['count_file'], row['qc_file'])
            
            counts['relative_counts'] = relative_counts
            mito_counts = counts[counts['Chr'].str.startswith('M94286')]
            
            if not mito_counts.empty:
                total = mito_counts['relative_counts'].sum()
                percentages.append(total)
        
        return np.mean(percentages) if percentages else 0
    
    def run_pipeline(self) -> None:
        """Execute the complete analysis pipeline."""
        logger.info("Starting VSG meta-analysis pipeline...")
        
        # Load data
        self.load_experiment_table()
        self.load_vsg_dictionary()
        self.validate_files()
        
        # Process samples
        self.identify_top_vsgs()
        
        # Generate configurations
        self.generate_experiment_config()
        
        # Save experiment configuration
        config_path = self.output_dir / 'exp_config.json'
        with open(config_path, 'w') as f:
            json.dump(self.exp_config, f, indent=4)
        logger.info(f"Saved experiment configuration to {config_path}")
        
        # Analyze silent VSGs
        silent_results = self.analyze_silent_vsgs()
        silent_path = self.output_dir / 'silentCsvData.csv'
        silent_results.to_csv(silent_path)
        logger.info(f"Saved silent VSG data to {silent_path}")
        
        # Analyze main VSGs
        main_results = self.analyze_main_vsgs()
        main_path = self.output_dir / 'mainCsvData.csv'
        main_results.to_csv(main_path)
        logger.info(f"Saved main VSG data to {main_path}")
        
        # Calculate QC metrics
        qc_results = self.calculate_qc_metrics()
        qc_path = self.output_dir / 'QC.csv'
        qc_results.to_csv(qc_path, index=False)
        logger.info(f"Saved QC metrics to {qc_path}")
        
        logger.info("Pipeline completed successfully!")
        
        # Print summary statistics
        self._print_summary(silent_results, main_results, qc_results)
    
    def _print_summary(self, silent_df: pd.DataFrame, main_df: pd.DataFrame, qc_df: pd.DataFrame) -> None:
        """Print summary statistics."""
        print("\n" + "="*50)
        print("PIPELINE SUMMARY")
        print("="*50)
        print(f"Total experiments processed: {len(self.exp_config)}")
        print(f"Total samples analyzed: {len(self.df)}")
        print(f"\nSilent VSG families detected: {silent_df.shape[1] - 1}")  # -1 for sum column
        print(f"Top 5 derepressed experiments (silent VSGs):")
        print(silent_df.nlargest(5, 'sum')[['sum']])
        print(f"\nTop 5 experiments with main VSG downregulation:")
        print(main_df.nsmallest(5, 'counts')[['counts']])
        print(f"\nQC Summary:")
        print(f"  Mean main VSG percentage: {qc_df['Main_VSG_perc'].mean():.2%}")
        print(f"  Mean mitochondrial percentage: {qc_df['Mito_perc'].mean():.2%}")
        print("="*50)


def main():
    """Main execution function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='🧬 VSG Meta-Analysis Pipeline - Analyze VSG expression across RNA-seq experiments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📊 Output Files:
  • exp_config.json    - Experiment configuration
  • silentCsvData.csv  - Silent VSG expression changes
  • mainCsvData.csv    - Main VSG expression changes  
  • QC.csv            - Quality control metrics

🌐 Visualization: https://vsgs-web-server.pages.dev/
📚 Documentation: https://github.com/mtinti/myRna-seq
        """
    )
    parser.add_argument(
        '--experiment-table',
        default='experiment_table.txt',
        help='Path to experiment metadata table'
    )
    parser.add_argument(
        '--vsg-dict',
        default='vsg_dic.txt',
        help='Path to VSG dictionary file'
    )
    parser.add_argument(
        '--base-path',
        default='myRna-seq/results/result_vsgs',
        help='Base path for result files'
    )
    parser.add_argument(
        '--output-dir',
        default='output',
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    # Print banner using ANSI codes directly
    print("\033[96m")  # Cyan color
    print("╔════════════════════════════════════════════════════════╗")
    print("║          🧬 VSG META-ANALYSIS PIPELINE 🧬              ║")
    print("║                                                        ║")
    print("║  Analyzing Variant Surface Glycoprotein Expression    ║")
    print("║           in Trypanosoma brucei                       ║")
    print("╚════════════════════════════════════════════════════════╝")
    print("\033[0m\n")  # Reset color
    
    # Initialize and run pipeline
    pipeline = VSGAnalysisPipeline(
        experiment_table_path=args.experiment_table,
        vsg_dict_path=args.vsg_dict,
        base_result_path=args.base_path,
        output_dir=args.output_dir
    )
    
    try:
        pipeline.run_pipeline()
    except FileNotFoundError as e:
        Console.error(f"Required file not found: {e}")
        Console.info("Please check that all input files exist and paths are correct")
        raise
    except Exception as e:
        Console.error(f"Pipeline failed: {e}")
        raise


if __name__ == '__main__':
    main()