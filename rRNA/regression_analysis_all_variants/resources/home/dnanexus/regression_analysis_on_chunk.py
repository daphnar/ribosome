import pandas as pd
import sys
import numpy as np
from scipy.stats import spearmanr,mannwhitneyu
import statsmodels.api as sm
import warnings
from statsmodels.tools.sm_exceptions import ConvergenceWarning

# coding_path="/home/dnanexus/mnt/rDNA Variations/coding19.tsv"
# field="/home/dnanexus/mnt/rDNA Variations/field.txt"
# covariates_path="/home/dnanexus/mnt/rDNA Variations/covariates.extra.csv"
# variants_path="/home/dnanexus/mnt/rDNA Variations/merged_rdna_variant_frequencies.unrelated_WB.csv"
# disease_path="/home/dnanexus/mnt/rDNA Variations/icd10_and_cancers.all.csv"
# # Read the CSV file and process the required chunk


def one_hot_encode_columns(df, categorical_columns):
    """
    Replaces categorical columns in the DataFrame with their one-hot encoded counterparts.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - categorical_columns (list of str): List of column names to be one-hot encoded.

    Returns:
    - pd.DataFrame: A new DataFrame with one-hot encoded columns replacing the original columns.
    """
    # Ensure the input columns are in the DataFrame
    for col in categorical_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame.")

    # Perform one-hot encoding
    df_encoded = pd.get_dummies(df, columns=categorical_columns, drop_first=True,dtype=float)
    return df_encoded


def get_covariates_and_variants():
    merged_rdna_variant_frequencies_path = 'merged_rdna_variant_frequencies.unrelated_WB.csv'
    merged_rdna_count_frequencies_path = 'merged_rdna_count_frequencies.all.csv'

    conversion = pd.read_csv('field.txt', sep='\t')
    covariates_df = pd.read_csv('covariates.extra.csv', index_col=0).set_index('eid')
    covariates_df.columns = [
        conversion[conversion['field_id'] == int(covariate.split('_')[0][1:])]['title'].values[0] + covariate[-3:] if
        covariate[0] == 'p' else covariate for covariate in covariates_df.columns]
    covariates_df = one_hot_encode_columns(covariates_df,
                                           ['Sexp31', 'Genotype measurement batch000', 'Release tranche050',
                                            'Shipment batch number053'])

    print('Got all covariates')
    all_variant_frequecies = pd.read_csv(merged_rdna_variant_frequencies_path).set_index(
        ['gene', 'position', 'nuc', 'variant']).dropna(how='all', axis=1)

    sample_names = list(all_variant_frequecies.columns.values)
    print('The sample size is %s, number of variants %s' % (len(sample_names), all_variant_frequecies.shape[0]))

    all_rdna_fraction = pd.read_csv(merged_rdna_count_frequencies_path, index_col=0)['Frac_47s'].to_frame()

    sample_names = list(set(all_rdna_fraction.index).intersection(sample_names))
    all_rdna_fraction = all_rdna_fraction.loc[sample_names]

    print('Starting to work on phenotype associations with %s UKBB samples' % len(sample_names))
    all_variant_frequecies = all_variant_frequecies.reset_index()
    all_variant_frequecies = all_variant_frequecies[all_variant_frequecies['variant'] != 'Reference'].set_index(
        ['gene', 'position', 'nuc', 'variant'])

    return covariates_df,all_variant_frequecies,all_rdna_fraction

def run_linear_regression(start_idx,end_idx,output_f):
    regression_output_path = output_f
    covariates_df,all_variant_frequecies,all_rdna_fraction=get_covariates_and_variants()
    variants=all_variant_frequecies.index
    phenotypes_df = pd.read_csv('phenotypes.csv', low_memory=False).set_index('eid')
    phenotypes_df = phenotypes_df.loc[:, phenotypes_df.columns.str.endswith('_i0')]  # MAYBE I AM REMOVING TOO MUCH
    phenotypes_df = pd.merge(phenotypes_df.reset_index(), covariates_df, on=['eid'], how='outer')

    print('Got all merged phenotypes')

    phenotypes_df['Sample']=phenotypes_df['eid'].apply(lambda x: str(x)+'_24048_0_0')
    phenotypes_df=phenotypes_df.set_index('Sample')

    print('Done renaming samples to match the rDNA naming')

    regression_output=[]
    phenotypes_df = phenotypes_df.join(all_rdna_fraction, how='outer')
    conversion = pd.read_csv('field.txt', sep='\t')
    phenotypes = phenotypes_df.columns.values[start_idx:end_idx + 1]
    sample_names = list(all_variant_frequecies.columns.values)

    for phenotype in phenotypes:
        try:
            phenotype_name=conversion[conversion['field_id']==int(phenotype[1:-3])]['title'].values[0]
        except:
            print('%s was not successfully converted'%phenotype)
            if phenotype=='Age' or phenotype=='Sex':
                phenotype_name=phenotype
            else:
                continue
        print('Working on %s'%phenotype_name)
        current_samples = phenotypes_df.loc[sample_names, phenotype].dropna().index
        df = phenotypes_df.loc[current_samples,list([phenotype,'Frac_47s'])+list(covariates_df.columns)].dropna()
        current_samples = df.index
        print('Number of samples for %s: %s'%(phenotype_name,len(current_samples)))
        X = df.drop(phenotype,axis=1)#[['Frac_47s']]
        X = sm.add_constant(X)
        y = df[phenotype]
        print('Model ready for copy number regression')
        try:
            model = sm.OLS(y, X).fit()
            print('Ran model')
            r = spearmanr(phenotypes_df.loc[current_samples, 'Frac_47s'], phenotypes_df.loc[current_samples, phenotype])[0]
            print('Ran correlation')
            regression_output.append([phenotype_name, 'copy number','','','',r]+list(model.pvalues.values)+['']+list(model.params.values)+[''])
        except ValueError:
            print('%s failed ValueError' % phenotype_name)
            pass
        except TypeError:
            print('%s failed TypeError' % phenotype_name)
            pass
        for idx,variant in enumerate(variants):
            r=spearmanr(all_variant_frequecies.loc[variant,current_samples],phenotypes_df.loc[current_samples,phenotype])[0]
            df.loc[current_samples,'rdna']=all_variant_frequecies.loc[variant,current_samples]
            X = df.drop(phenotype,axis=1)#
            X = sm.add_constant(X)
            y = df[phenotype]
            model = sm.OLS(y, X).fit()
            print('model done %s'%str(variant))
            gene=variant[0]
            position=variant[1]
            nuc=variant[2]
            variant_str=variant[3]
            regression_output.append([phenotype_name, gene, position, nuc, variant_str, r] + list(model.pvalues.values) + list(model.params.values))

        print('Wrote stats DF', regression_output_path)
        pd.DataFrame(regression_output,columns=['phenotype','gene','position','nuc','variant','correlation']+[name+'_pvalue' for name in model.params.index]+list(model.pvalues.index)).to_csv(regression_output_path)


def run_logistic_regression(start_idx,end_idx,output_f,control_type='any'):
    covariates_df,all_variant_frequecies,all_rdna_fraction=get_covariates_and_variants()

    phenotypes_df = pd.read_csv('icd10_and_cancers.all.csv', low_memory=False).set_index('eid')

    diseases=phenotypes_df['p41270']
    disease_list=[]
    for sample in diseases.index:
        if type(diseases[sample])==str:
            disease_list+=eval(diseases[sample])
    disease_list=pd.Series(disease_list).value_counts()
    disease_list=disease_list[disease_list>1000]
    disease_list=disease_list[start_idx:end_idx+1]
    phenotypes_df = phenotypes_df.loc[:, 'p41270'].fillna('[]').apply(
        lambda x: eval(x) if x else []).to_frame()  # Cancer field
    phenotypes_df = pd.merge(phenotypes_df.reset_index(), covariates_df, on=['eid'], how='outer').dropna()
    phenotypes_df = phenotypes_df.reset_index()
    print('Got all merged phenotypes')

    phenotypes_df['Sample'] = phenotypes_df['eid'].apply(lambda x: str(x) + '_24048_0_0')
    phenotypes_df = phenotypes_df.set_index('Sample')

    variants = all_variant_frequecies.index.values

    conversion = pd.read_csv('coding19.tsv', sep='\t', index_col=0)
    phenotype_names = {k: conversion.loc[k, 'meaning'] for k in disease_list.index}

    regression_output_path = output_f
    regression_output = []
    phenotypes_df = phenotypes_df.join(all_rdna_fraction, how='outer')
    have_no_disease = phenotypes_df.loc[
        (phenotypes_df['p41270'].isna()) | (phenotypes_df['p41270'].apply(lambda x: x == []))].index

    warnings.simplefilter("error", ConvergenceWarning)
    for phenotype, phenotype_name in phenotype_names.items():  # phenotypes_df.columns:#'BMI','Height']:#['BMI','Height','Age','Sex','Frac_47s','Frac_coding','Frac_noncoding']:
        print('Working on %s' % phenotype_name)
        current_phenotypes = phenotypes_df.copy()
        disease_samples = current_phenotypes.dropna()
        disease_samples = disease_samples[disease_samples['p41270'].apply(lambda x: phenotype in x)].index
        if control_type=='any':
            current_phenotypes[phenotype] = 0
        else:
            current_phenotypes[phenotype] = np.nan
        current_phenotypes.loc[disease_samples, phenotype] = 1
        if control_type!='any':
            current_phenotypes.loc[current_phenotypes.index.isin(have_no_disease), phenotype] = 0
        print(len(disease_samples), ' samples with disease')
        current_phenotypes = current_phenotypes.loc[current_phenotypes[phenotype].dropna().index].dropna()
        current_samples = current_phenotypes.index
        print(current_phenotypes.columns)
        print(phenotype,covariates_df.columns)
        df = current_phenotypes.loc[current_samples, list([phenotype, 'Frac_47s']) + list(
            covariates_df.columns)].dropna()  ##this is needed to get rid of unwanted columns
        current_samples = df.index
        print('Number of samples for %s: %s' % (phenotype_name, len(current_samples)))
        X = df.drop(phenotype, axis=1).copy()
        X = (X - X.mean()) / X.std()
        X = sm.add_constant(X)
        y = df[phenotype]  # y has to be 0,1
        print('Model ready for copy number regression')
        r = spearmanr(df.loc[current_samples, 'Frac_47s'], df.loc[current_samples, phenotype])[0]
        try:
            model = sm.Logit(y, X).fit()
            print('Ran model')
            print('Ran correlation')
            regression_output.append(
                [phenotype_name, 'copy number', '', '', '', r, ''] + list(model.pvalues.values) + [''] + list(
                    model.params.values) + [''])
        except Exception:
            print('%s did not converge!' % phenotype_name)
            #regression_output.append([phenotype_name,  'copy number', '', '', '', r, ''] + ['@'] * 324)#324=columns

        for idx, variant in enumerate(variants):
            r = spearmanr(all_variant_frequecies.loc[variant, current_samples], df.loc[current_samples, phenotype])[0]
            df_current = df.dropna().copy()
            df_current.loc[:, 'rdna'] = all_variant_frequecies.loc[variant, df_current.index]
            ranksums_p = mannwhitneyu(df_current.loc[df_current[phenotype]==1,'rdna'],
                                      df_current.loc[df_current[phenotype]==0,'rdna'])[1]
            X = df_current.drop(phenotype, axis=1)
            X = (X - X.mean()) / X.std()
            X = sm.add_constant(X)
            y = df_current[phenotype]
            gene = variant[0]
            position = variant[1]
            nuc = variant[2]
            variant_str = variant[3]
            try:
                model = sm.Logit(y, X).fit()
                print('model done %s' % str(variant))
                regression_output.append([phenotype_name, gene, position, nuc, variant_str, r, ranksums_p] + list(
                    model.pvalues.values) + list(model.params.values))
            except Exception:
                print('%s did not converge!' % str(variant))
                #regression_output.append([phenotype_name, gene, position, nuc, variant_str, r, ranksums_p] + ['@']*324)
                pass
        print('Wrote stats DF', regression_output_path)
        try:
            pd.DataFrame(regression_output,
                     columns=['phenotype', 'gene', 'position', 'nuc', 'variant', 'correlation', 'ranksums_p'] + [
                         name + '_pvalue' for name in model.params.index] + list(model.pvalues.index)).to_csv(
            regression_output_path)
        except Exception:
            print('%s model failed for all!!!!!' % phenotype)
            pass

if __name__=='__main__':
    regression_type = sys.argv[1]
    start_idx = int(sys.argv[2])
    end_idx = int(sys.argv[3])
    output_f= sys.argv[4]
    control_type = sys.argv[5]
    if regression_type=='logistic':
        run_logistic_regression(start_idx,end_idx,output_f,control_type)
    elif regression_type=='linear':
        run_linear_regression(start_idx,end_idx,output_f)

