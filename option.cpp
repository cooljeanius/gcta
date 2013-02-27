/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * GCTA options
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void option(int option_num, char* option_str[]);

int main(int argc, char* argv[])
{
  //  argc=8;
   // char * x[8]={"gcta", "--bfile", "K:/scratch/JointMA/GIANT2/GIANT_180_HT_SNP", "--massoc-file", "K:/scratch/JointMA/GIANT2/GIANT_HT_3M_SNPs_QCed.ma", "--chr", "11", "--massoc-joint"};
  //  argv=x;
 /*   string tmp="1";
    vector<string> str_buf;
    cout<<StrFunc::split_str(tmp, str_buf)<<" strings"<<endl;
    int m=0;
    for( m=0; m<str_buf.size(); m++)cout<<str_buf[m]<<endl;
    cin>>m;
    exit(0);
*/

	cout<<"*******************************************************************"<<endl;
	cout<<"* Genome-wide Complex Trait Analysis (GCTA)"<<endl;
	cout<<"* version 0.93.9"<<endl;
	cout<<"* (C) 2010 Jian Yang, Hong Lee, Michael Goddard and Peter Visscher"<<endl;
	cout<<"* GNU General Public License, v2"<<endl;
	cout<<"* Queensland Institute of Medical Research"<<endl;
	cout<<"*******************************************************************"<<endl;
	long int time_used=0, start=time(NULL);
	time_t curr=time(0);
	cout<<"Analysis started: "<<ctime(&curr)<<endl;
	cout<<"Options:"<<endl;
	try{ option(argc, argv); }
	catch(const string &err_msg){ cerr<<err_msg<<endl; }
	catch(const char *err_msg){ cerr<<err_msg<<endl; }
	curr=time(0);
	cout<<"\nAnalysis finished: "<<ctime(&curr);
	time_used=time(NULL)-start;
	cout<<"Computational time: "<<time_used/3600<<":"<<(time_used%3600)/60<<":"<<time_used%60<<endl;

	return 0;
}

void option(int option_num, char* option_str[])
{
	int i=0, j=0;

	// initialize parameters for the genome() function
/*	bool genome_flag=false;
	string genome_popSize("10000"), genome_rec("0.0001");
	int genome_numPieces=100, genome_pieceLen=10000, genome_numIndepRegion=1, genome_SNP=-1;
	int genome_nSubPOP=2;
	vector<int> genome_nSubSample(genome_nSubPOP);
	for(i=0; i<genome_nSubPOP; i++) genome_nSubSample[i]=10;
	double genome_mut=1e-8, genome_mig=2.5e-4;
*/
    // raw genotype data
    string RG_fname_file="", RG_summary_file="";
    double GC_cutoff=0.7;

	// data management
	string bfile="", bfile2="", update_sex_file="", update_freq_file="", update_refA_file="", kp_indi_file="", rm_indi_file="", extract_snp_file="", exclude_snp_file="", extract_snp_name="", exclude_snp_name="", out="gcta";
	bool SNP_major=false, bfile_flag=false, make_bed_flag=false, dose_mach_flag=false, dose_beagle_flag=false, bfile2_flag=false, out_freq_flag=false, out_ssq_flag=false;
	bool ref_A=false, recode=false, recode_nomiss=false, save_ram=false, autosome_flag=false;
	int autosome_num=22, extract_chr_start=0, extract_chr_end=0;
	string dose_file="", dose_info_file="", update_dose_rsq_file="";
	double maf=0.0, max_maf=0.0, dose_Rsq_cutoff=0.0;

	// GRM
	bool ibc=false, ibc_all=false, grm_flag=false, m_grm_flag=false, make_grm_flag=false, make_grm_xchar_flag=false, grm_out_bin_flag=false;
	bool pca_flag=false;
	double grm_adj_fac=-2.0, grm_cutoff=-2.0, rm_high_ld_cutoff=-1.0;
	int dosage_compen=-2, out_pc_num=20, grm_mtd=0;
	string grm_file="", paa_file="";

	// initialize paramters for estimating LD structure
	string LD_file="";
	bool LD=false, LD_search=false, LD_i=false;
	int LD_step=10;
	double LD_window=1e7, LD_sig=0.05;

	// initialize paramters for simulation based on real genotype data
	bool simu_qt_flag=false, simu_cc=false, simu_emb_flag=false, simu_output_causal=false;
	int simu_rep=1, simu_case_num=0, simu_control_num=0;
	double simu_h2=0.1, simu_K=0.1, simu_gener=100;
	string simu_causal="";

	// estimate genetic distance based on hapmap_data
	bool hapmap_genet_dst=false;
	string hapmap_genet_dst_file="";

	// REML analysis
	int mphen=1, reml_mtd=0, MaxIter=100;
	double prevalence=-2.0;
	bool reml_flag=false, pred_rand_eff=false, est_fix_eff=false, blup_snp_flag=false, no_constrain=false, reml_lrt_flag=false, no_lrt=false;
	string phen_file="", qcovar_file="", covar_file="", qgxe_file="", gxe_file="", blup_indi_file="";
	vector<double> reml_priors, reml_priors_var;
	vector<int> reml_drop;
	reml_drop.push_back(1);

    int argc=option_num;
    vector<char *> argv(option_num+2);
    for(i=0; i<option_num; i++) argv[i]=option_str[i];
    argv[option_num]="gcta"; argv[option_num+1]="gcta";
	for(i=1; i<argc; i++){
        // raw genotype data
		if(strcmp(argv[i],"--raw-files")==0){
			RG_fname_file=argv[++i];
			cout<<"--raw-files "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--raw-summary")==0){
			RG_summary_file=argv[++i];
			cout<<"--raw-summary "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--gencall")==0){
			GC_cutoff=atof(argv[++i]);
			cout<<"--gencall "<<GC_cutoff<<endl;
			if(GC_cutoff<0.0 || GC_cutoff>1.0) throw("\nError: --gencall should be within the range from 0 to 1.\n");
		}

		// data management
		else if(strcmp(argv[i],"--bfile")==0){
			bfile_flag=true;
			bfile=argv[++i];
			cout<<"--bfile "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--make-bed")==0){
			make_bed_flag=true;
			cout<<"--make-bed "<<endl;
		}
		else if(strcmp(argv[i],"--bfile2")==0){
			bfile2_flag=true;
			bfile2=argv[++i];
			cout<<"--bfile2 "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--dosage-mach")==0){
			dose_mach_flag=true;
			dose_beagle_flag=false;
			dose_file=argv[++i];
			dose_info_file=argv[++i];
			out=dose_file;
			cout<<"--dosage-mach "<<dose_file<<" "<<dose_info_file<<endl;
		}
		else if(strcmp(argv[i],"--dosage-beagle")==0){
			dose_beagle_flag=true;
			dose_mach_flag=false;
			dose_file=argv[++i];
			dose_info_file=argv[++i];
			out=dose_file;
			cout<<"--dosage-beagle "<<dose_file<<" "<<dose_info_file<<endl;
		}
		else if(strcmp(argv[i],"--imput-rsq")==0){
			dose_Rsq_cutoff=atof(argv[++i]);
			cout<<"--imput-rsq "<<dose_Rsq_cutoff<<endl;
			if(dose_Rsq_cutoff<0.0 || dose_Rsq_cutoff>1.0) throw("\nError: --imput-rsq should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--update-imput-rsq")==0){
			update_dose_rsq_file=argv[++i];
			cout<<"--update-imput-rsq "<<update_dose_rsq_file<<endl;
		}
		else if(strcmp(argv[i],"--update-freq")==0){
			update_freq_file=argv[++i];
			cout<<"--update-freq "<<update_freq_file<<endl;
		}
		else if(strcmp(argv[i],"--update-ref-allele")==0){
			update_refA_file=argv[++i];
			cout<<"--update-ref-allele "<<update_refA_file<<endl;
		}
		else if(strcmp(argv[i],"--keep")==0){
			kp_indi_file=argv[++i];
			cout<<"--keep "<<kp_indi_file<<endl;
		}
		else if(strcmp(argv[i],"--remove")==0){
			rm_indi_file=argv[++i];
			cout<<"--remove "<<rm_indi_file<<endl;
		}
		else if(strcmp(argv[i],"--update-sex")==0){
			update_sex_file=argv[++i];
			cout<<"--update-sex "<<update_sex_file<<endl;
		}
		else if(strcmp(argv[i],"--chr")==0){
			extract_chr_start=extract_chr_end=atoi(argv[++i]);
			cout<<"--chr "<<extract_chr_start<<endl;
			if(extract_chr_start<1 || extract_chr_start>100) throw("\nError: --chr should be within the range from 1 to 100.\n");
		}
		else if(strcmp(argv[i],"--autosome-num")==0){
			autosome_num=atoi(argv[++i]);
			cout<<"--autosome-num "<<autosome_num<<endl;
            if(autosome_num<1 || autosome_num>100) throw("\nError: invalid number specified after the option --autosome-num.\n");
		}
		else if(strcmp(argv[i],"--autosome")==0){
		    autosome_flag=true;
			cout<<"--autosome"<<endl;
		}
		else if(strcmp(argv[i],"--extract")==0){
			extract_snp_file=argv[++i];
			cout<<"--extract "<<extract_snp_file<<endl;
		}
		else if(strcmp(argv[i],"--exclude")==0){
			exclude_snp_file=argv[++i];
			cout<<"--exclude "<<exclude_snp_file<<endl;
		}
		else if(strcmp(argv[i],"--extract-snp")==0){
			extract_snp_name=argv[++i];
			cout<<"--extract-snp "<<extract_snp_name<<endl;
		}
		else if(strcmp(argv[i],"--exclude-snp")==0){
			exclude_snp_name=argv[++i];
			cout<<"--exclude-snp "<<exclude_snp_name<<endl;
		}		
        else if(strcmp(argv[i],"--maf")==0){
			maf=atof(argv[++i]);
			cout<<"--maf "<<maf<<endl;
			if(maf<0 || maf>0.5) throw("\nError: --maf should be within the range from 0 to 0.5.\n");
		}
		else if(strcmp(argv[i],"--max-maf")==0){
			max_maf=atof(argv[++i]);
			cout<<"--max-maf "<<max_maf<<endl;
			if(max_maf<0 || max_maf>0.5) throw("\nError: --max-maf should be within the range from 0 to 0.5.\n");
		}
		else if(strcmp(argv[i],"--out")==0){
			out=argv[++i];
			cout<<"--out "<<out<<endl;
		}
		else if(strcmp(argv[i],"--freq")==0){
			out_freq_flag=true;
			cout<<"--freq"<<endl;
		}
		else if(strcmp(argv[i],"--ssq")==0){
			out_ssq_flag=true;
			cout<<"--ssq"<<endl;
		}
		else if(strcmp(argv[i],"--recode")==0){
			recode=true;
			cout<<"--recode"<<endl;
		}
		else if(strcmp(argv[i],"--recode-nomiss")==0){
			recode_nomiss=true;
			cout<<"--recode-nomiss"<<endl;
		}
		else if(strcmp(argv[i],"--save-ram")==0){
			save_ram=true;
			cout<<"--save-ram"<<endl;
		}
		// GRM
		else if(strcmp(argv[i],"--paa")==0){
			paa_file=argv[++i];
			cout<<"--paa "<<paa_file<<endl;
		}
		else if(strcmp(argv[i],"--ibc")==0){
			ibc=true;
			cout<<"--ibc"<<endl;
		}
		else if(strcmp(argv[i],"--ibc-all")==0){
		    ibc=ibc_all=true;
			cout<<"--ibc-all"<<endl;
		}
		else if(strcmp(argv[i],"--mgrm")==0){
			m_grm_flag=true;
			grm_file=argv[++i];
			cout<<"--mgrm "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--grm")==0){
			grm_flag=true;
			grm_file=argv[++i];
			cout<<"--grm "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--grm-alg")==0){
			grm_mtd=atoi(argv[++i]);
			cout<<"--grm-alg "<<grm_mtd<<endl;
			if(grm_mtd<0 || grm_mtd>1) throw("\nError: --grm-alg should be 0 or 1.\n");
		}
		else if(strcmp(argv[i],"--rm-high-ld")==0){
			rm_high_ld_cutoff=atof(argv[++i]);
			cout<<"--rm-high-ld "<<rm_high_ld_cutoff<<endl;
			if(rm_high_ld_cutoff<=0 || rm_high_ld_cutoff>=1) throw("\nError: the value to be specified after --rm-high-ld should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--make-grm")==0){
			make_grm_flag=true;
			cout<<"--make-grm"<<endl;
		}
		else if(strcmp(argv[i],"--make-grm-xchr")==0){
			make_grm_xchar_flag=true;
			cout<<"--make-grm-xchr"<<endl;
		}
		else if(strcmp(argv[i],"--make-grm-bin")==0){
		    make_grm_flag=true;
			grm_out_bin_flag=true;
			cout<<"--make-grm-bin"<<endl;
		}
		else if(strcmp(argv[i],"--make-grm-xchr-bin")==0){
			make_grm_xchar_flag=true;
			grm_out_bin_flag=true;
			cout<<"--make-grm-xchr-bin"<<endl;
		}
		else if(strcmp(argv[i],"--grm-adj")==0){
			grm_adj_fac=atof(argv[++i]);
			cout<<"--grm-adj "<<grm_adj_fac<<endl;
			if(grm_adj_fac<0 || grm_adj_fac>1) throw("\nError: the value to be specified after --grm-adj should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--dc")==0){
			dosage_compen=atoi(argv[++i]);
			cout<<"--dc "<<dosage_compen<<endl;
			if(dosage_compen!=0 && dosage_compen!=1) throw("\nError: the value to be specified after --dc should be 0 or 1.\n");
		}
		else if(strcmp(argv[i],"--grm-cutoff")==0){
			grm_cutoff=atof(argv[++i]);
			cout<<"--grm-cutoff "<<grm_cutoff<<endl;
			if(grm_cutoff<0 || grm_cutoff>1) throw("\nError: the value to be specified after --grm-cutoff should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--pca")==0){
		    pca_flag=true;
            i++;
            if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) { out_pc_num=20; i--; }
            else out_pc_num=atoi(argv[i]);
			cout<<"--pca "<<out_pc_num<<endl;
			if(out_pc_num<1) throw("\nError: the value to be specified after --pca should be positive.\n");
		}
		// estimation of LD structure
		else if(strcmp(argv[i],"--ld")==0){
			LD=true;
			LD_file=argv[++i];
			cout<<"--ld "<<LD_file<<endl;
		}
		else if(strcmp(argv[i],"--ld-step")==0){
			LD_search=true;
			LD_step=atoi(argv[++i]);
			cout<<"--ld-step "<<LD_step<<endl;
			if(LD_step<1 || LD_step>20) throw("\nError: --ld-step should be within the range from 1 to 20.\n");
		}
		else if(strcmp(argv[i],"--ld-wind")==0){
			LD_window=atof(argv[++i]);
			cout<<"--ld-wind "<<LD_window<<endl;
			LD_window*=1000;
			if(LD_window<1e3 || LD_window>2e8) throw("\nError: --ld-wind should be within the range from 1Kb to 20Mb.\n");
		}
		else if(strcmp(argv[i],"--ld-sig")==0){
			LD_sig=atof(argv[++i]);
			cout<<"--ld-sig "<<LD_sig<<endl;
			if(LD_sig<=0) throw("\nError: --ld-sig should be > 0.\n");
		}
		else if(strcmp(argv[i],"--ld-i")==0){
			LD_i=true;
			cout<<"--ld-i"<<endl;
		}
		// simulation based on real genotype data
		else if(strcmp(argv[i],"--simu-qt")==0){
			simu_qt_flag=true;
			cout<<"--simu-qt"<<endl;
		}
		else if(strcmp(argv[i],"--simu-cc")==0){
			simu_cc=true;
			simu_case_num=atoi(argv[++i]);
			simu_control_num=atoi(argv[++i]);
			cout<<"--simu-cc "<<simu_case_num<<" "<<simu_control_num<<endl;
            if(simu_case_num<10) throw("Error: --simu-cc, Invalid number of cases. Minimun number 10.");
            if(simu_control_num<10) throw("Error: --simu-cc, Invalid number of controls. Minimum number 10.");
		}
		else if(strcmp(argv[i],"--simu-rep")==0){
			simu_rep=atoi(argv[++i]);
			cout<<"--simu-rep "<<simu_rep<<endl;
            if(simu_rep<1 || simu_rep>10000) throw("Error: --simu-rep should be within the range from 1 to 10000.");
		}
		else if(strcmp(argv[i],"--simu-hsq")==0){
			simu_h2=atof(argv[++i]);
			cout<<"--simu-hsq "<<simu_h2<<endl;
            if(simu_h2>1.0 || simu_h2<0.0) throw("Error: --simu-h2 should be within the range from 0 to 1.");
		}
		else if(strcmp(argv[i],"--simu-k")==0){
			simu_K=atof(argv[++i]);
			cout<<"--simu-k "<<simu_K<<endl;
            if(simu_K>0.5 || simu_K<0.0001) throw("Error: --simu-K should be within the range from 0.0001 to 0.5.");
		}
		else if(strcmp(argv[i],"--simu-causal-loci")==0){
			simu_causal=argv[++i];
			cout<<"--simu-causal-loci "<<simu_causal<<endl;
		}
		else if(strcmp(argv[i],"--simu-embayesb")==0){ // internal
			simu_emb_flag=true;
			cout<<"--simu-embayesb"<<endl;
		}
		else if(strcmp(argv[i],"--simu-ouput-causal")==0){ // internal
			simu_output_causal=true;
			cout<<"--simu-output-causal"<<endl;
		}
/*		else if(strcmp(argv[i],"--simu-gener")==0){
			simu_gener=atoi(argv[++i]);
            if(simu_gener<0 || simu_gener>1e5) throw("Error: --simu-gener should be within the range from 0 to 100000.");
			cout<<"--simu-gener "<<simu_gener<<endl;
		}
		*/
		// calculate genetic dst based on HapMap data
		else if(strcmp(argv[i],"--hapmap-genet-dst")==0){
			hapmap_genet_dst=true;
			hapmap_genet_dst_file=argv[++i];
			cout<<"--hapmap-genet-dst "<<hapmap_genet_dst_file<<endl;
		}
		// estimate variance explained by all SNPs
		else if(strcmp(argv[i],"--reml")==0){
			reml_flag=true;
			cout<<"--reml"<<endl;
		}
		else if(strcmp(argv[i],"--prevalence")==0){
			prevalence=atof(argv[++i]);
			cout<<"--prevalence "<<prevalence<<endl;
			if(prevalence<=0 || prevalence>=1) throw("\nError: --prevalence should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--reml-pred-rand")==0){
			pred_rand_eff=true;
			cout<<"--reml-pred-rand"<<endl;
		}
		else if(strcmp(argv[i],"--reml-est-fix")==0){
			est_fix_eff=true;
			cout<<"--reml-est-fix"<<endl;
		}
		else if(strcmp(argv[i],"--reml-alg")==0){
			reml_mtd=atoi(argv[++i]);
			cout<<"--reml-alg "<<reml_mtd<<endl;
			if(reml_mtd<0 || reml_mtd>2) throw("\nError: --reml-alg should be 0, 1 or 2.\n");
		}
		else if(strcmp(argv[i],"--reml-no-constrain")==0){
			no_constrain=true;
			cout<<"--reml-no-constrain"<<endl;
		}
		else if(strcmp(argv[i],"--reml-priors")==0){
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_priors.push_back(atof(argv[i]));
		    }
		    i--;
			cout<<"--reml-priors ";
			bool err_flag=false;
			for(j=0; j<reml_priors.size(); j++){
			    cout<<reml_priors[j]<<" ";
			    if(reml_priors[j]>1.0 || reml_priors[j]<0.0) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_priors.empty()) throw("\nError: --reml-priors. Prior values should be within the range from 0 to 1.\n");
			if(CommFunc::sum(reml_priors)>1.0) throw("\nError: --reml-priors. The sum of all prior values should not exceed 1.0.");
		}
		else if(strcmp(argv[i],"--reml-priors-var")==0){
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_priors_var.push_back(atof(argv[i]));
		    }
		    i--;
			cout<<"--reml-priors-var ";
			bool err_flag=false;
			for(j=0; j<reml_priors_var.size(); j++){
			    cout<<reml_priors_var[j]<<" ";
			    if(reml_priors_var[j]<0.0) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_priors_var.empty()) throw("\nError: --reml-priors-var. Prior values should be positive.\n");
		}
		else if(strcmp(argv[i],"--reml-no-lrt")==0){
			no_lrt=true;
			cout<<"--reml-no-lrt"<<endl;
		}
		else if(strcmp(argv[i],"--reml-lrt")==0){
            reml_lrt_flag=true;
		    reml_drop.clear();
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_drop.push_back(atoi(argv[i]));
		    }
		    i--;
			cout<<"--reml-lrt ";
			bool err_flag=false;
			for(j=0; j<reml_drop.size(); j++){
			    cout<<reml_drop[j]<<" ";
			    if(reml_drop[j]<1) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_drop.empty()) throw("\nError: invalid values specified after --reml-lrt.\n");
		}
		else if(strcmp(argv[i],"--reml-maxit")==0){
			MaxIter=atoi(argv[++i]);
			cout<<"--reml-maxit "<<MaxIter<<endl;
			if(MaxIter<1 || MaxIter>10000) throw("\nError: --reml-maxit should be within the range from 1 to 10000.\n");
		}
		else if(strcmp(argv[i],"--pheno")==0){
			phen_file=argv[++i];
			cout<<"--pheno "<<phen_file<<endl;
		}
		else if(strcmp(argv[i],"--mpheno")==0){
			mphen=atoi(argv[++i]);
			cout<<"--mpheno "<<mphen<<endl;
			if(mphen<1) throw("Error: --mpheno should be > 0.");
		}
		else if(strcmp(argv[i],"--qcovar")==0){
			qcovar_file=argv[++i];
			cout<<"--qcovar "<<qcovar_file<<endl;
		}
		else if(strcmp(argv[i],"--covar")==0){
			covar_file=argv[++i];
			cout<<"--covar "<<covar_file<<endl;
		}
		else if(strcmp(argv[i],"--gxqe")==0){
			qgxe_file=argv[++i];
			cout<<"--gxqe "<<qgxe_file<<endl;
		}
		else if(strcmp(argv[i],"--gxe")==0){
			gxe_file=argv[++i];
			cout<<"--gxe "<<gxe_file<<endl;
		}
		else if(strcmp(argv[i],"--blup-snp")==0){
		    blup_snp_flag=true;
			blup_indi_file=argv[++i];
			cout<<"--blup-snp "<<blup_indi_file<<endl;
		}
		else if(strcmp(argv[i],"gcta")==0) break;
		else{ stringstream errmsg; errmsg<<"\nError: invalid option \""<<argv[i]<<"\".\n"; throw(errmsg.str()); }
		// genome() function
/*		else if(strcmp(argv[i],"--genome-pop")==0){
			genome_flag=true;
			genome_nSubPOP=atoi(argv[i+1]);
			genome_nSubSample.clear();
			genome_nSubSample.resize(genome_nSubPOP);
			for(j=0; j<genome_nSubPOP; j++) genome_nSubSample[j]=atoi(argv[i+2+j]);
		}
		else if(strcmp(argv[i],"--genome-N")==0) genome_popSize.assign(argv[i+1]);
		else if(strcmp(argv[i],"--genome-c")==0) genome_numIndepRegion=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-pieces")==0) genome_numPieces=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-len")==0) genome_pieceLen=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-s")==0) genome_SNP=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-rec")==0) genome_rec.assign(argv[i+1]);
		else if(strcmp(argv[i],"--genome-mut")==0) genome_mut=atof(argv[i+1]);
		else if(strcmp(argv[i],"--genome-mig")==0) genome_mig=atof(argv[i+1]);
*/
	}
	// conflicted options
	cout<<endl;
	if(bfile2_flag && !bfile_flag) throw("Error: the option --bfile2 should always go with the option --bfile.");
	if(m_grm_flag){
	    if(grm_flag){ grm_flag=false; cout<<"Warning: --grm option suppressed by the --mgrm option."<<endl; }
	    if(grm_cutoff>-1.0){ grm_cutoff=-2.0; cout<<"Warning: --grm-cutoff option suppressed by the --mgrm option."<<endl; }
	}
	if(pca_flag){
	    if(grm_adj_fac>-1.0){ grm_adj_fac=-2.0; cout<<"Warning: --grm-adj option suppressed by the --pca option."<<endl; }
	    else if(dosage_compen>-1){ grm_adj_fac=-2; cout<<"Warning: --dosage-compen option suppressed by the --pca option."<<endl; }
	}
    if(!gxe_file.empty() && !grm_flag && !m_grm_flag){
        cout<<"Warning: --gxe option is ignored because there is no --grm or --mgrm option specified."<<endl;
        gxe_file="";
    }
    if(pred_rand_eff && !grm_flag && !m_grm_flag){
        cout<<"Warning: --reml-pred-rand option is ignored because there is no --grm or --mgrm option specified."<<endl;
        pred_rand_eff=false;
    }
    if(dosage_compen>-1 && update_sex_file.empty()) throw("Error: you need to specify the sex information for the individuals by the option --update-sex because of the option --dc.");
    if(bfile2_flag && update_freq_file.empty()) throw("Error: you need to update the allele frequency by the option --update-freq because there are two datasets.");

    // set autosome
    if(autosome_flag){
        extract_chr_start=1;
        extract_chr_end=autosome_num;
    }
    if(make_grm_xchar_flag) extract_chr_start=extract_chr_end=(autosome_num+1);

	// Implement
	cout<<endl;
    gcta *pter_gcta=new gcta(autosome_num, rm_high_ld_cutoff, out);//, *pter_gcta2=new gcta(autosome_num, rm_high_ld_cutoff, out);
    if(!RG_fname_file.empty()){
		if(RG_summary_file.empty()) throw("Error: please input the summary information for the raw data files by the option --raw-summary.");
        pter_gcta->read_IRG_fnames(RG_summary_file, RG_fname_file, GC_cutoff);
    }
    else if(bfile_flag){
		if(hapmap_genet_dst) pter_gcta->genet_dst(bfile, hapmap_genet_dst_file);
		else{
		    if(bfile2_flag){
		        cout<<"There are two datasets specified (in PLINK binary PED format).\nReading dataset 1 ..."<<endl;
		        if(update_freq_file.empty()) throw("Error: since there are two dataset, you should update the allele frequencies that are calculated in the combined dataset.");
		    }
			pter_gcta->read_famfile(bfile+".fam");
            if(!kp_indi_file.empty()) pter_gcta->keep_indi(kp_indi_file);
			if(!rm_indi_file.empty()) pter_gcta->remove_indi(rm_indi_file);
			if(!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
			if(!blup_indi_file.empty()) pter_gcta->read_indi_blup(blup_indi_file);
			pter_gcta->read_bimfile(bfile+".bim");
			if(!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
			if(!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
			if(extract_chr_start>0) pter_gcta->extract_chr(extract_chr_start, extract_chr_end);
			if(!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
			if(!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
			if(!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
			if(LD) pter_gcta->read_LD_target_SNPs(LD_file);
			pter_gcta->read_bedfile(bfile+".bed");
			if(!update_dose_rsq_file.empty()) pter_gcta->update_rsq(update_dose_rsq_file);
			if(!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
            if(dose_Rsq_cutoff>0.0) pter_gcta->filter_dosage_rsq(dose_Rsq_cutoff);
			if(maf>0) pter_gcta->filter_snp_maf(maf);
			if(max_maf>0.0) pter_gcta->filter_snp_max_maf(max_maf);
/*			if(bfile2_flag){
			    cout<<"Reading dataset 2 ..."<<endl;
			    pter_gcta2=new gcta();
                pter_gcta2->read_famfile(bfile2+".fam");
                if(!kp_indi_file.empty()) pter_gcta2->keep_indi(kp_indi_file);
                if(!rm_indi_file.empty()) pter_gcta2->remove_indi(rm_indi_file);
                if(!update_sex_file.empty()) pter_gcta2->update_sex(update_sex_file);
                pter_gcta2->read_bimfile(bfile2+".bim");
                pter_gcta2->read_bedfile(bfile2+".bed");
                pter_gcta->update_ref_A(update_refA_file);
			}*/
			if(out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
			else if(!paa_file.empty()) pter_gcta->paa(paa_file);
			else if(ibc) pter_gcta->ibc(ibc_all);
			else if(make_grm_flag) pter_gcta->make_grm(false, grm_out_bin_flag, grm_mtd);
			else if(make_grm_xchar_flag) pter_gcta->make_grm(true, grm_out_bin_flag, grm_mtd);
			else if(recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
			else if(LD) pter_gcta->LD_Blocks(LD_step, LD_window, LD_sig, LD_i, save_ram);
			else if(blup_snp_flag) pter_gcta->blup_snp_geno();
            else if(simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_output_causal, simu_emb_flag);
			else if(make_bed_flag) pter_gcta->save_plink();
		}
	}
	else if(dose_beagle_flag || dose_mach_flag){
		if(dose_mach_flag) pter_gcta->read_imp_info_mach(dose_info_file);
		else if(dose_beagle_flag) pter_gcta->read_imp_info_beagle(dose_info_file);
        if(!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
        if(!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
        if(!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
        if(!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
        if(extract_chr_start>0) cout<<"Warning: the option --chr, --autosome or --nonautosome is inactive for dosage data."<<endl;
        if(!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
        if(dose_mach_flag) pter_gcta->read_imp_dose_mach(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
		else if(dose_beagle_flag) pter_gcta->read_imp_dose_beagle(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
		if(!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
		if(!update_dose_rsq_file.empty()) pter_gcta->update_rsq(update_dose_rsq_file);
		if(!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
		if(dose_Rsq_cutoff>0.0) pter_gcta->filter_dosage_rsq(dose_Rsq_cutoff);
		if(maf>0.0) pter_gcta->filter_snp_maf(maf);
		if(max_maf>0.0) pter_gcta->filter_snp_max_maf(max_maf);
		if(out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
		if(make_grm_flag) pter_gcta->make_grm(false, grm_out_bin_flag, grm_mtd);
		else if(make_grm_xchar_flag) pter_gcta->make_grm(true, grm_out_bin_flag, grm_mtd);
		else if(recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
		else if(blup_snp_flag) pter_gcta->blup_snp_dosage();
        else if(simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_output_causal, simu_emb_flag);
	}
	else if(reml_flag){
		if(phen_file.empty()) throw("\nError: phenotype file is required for reml analysis.\n");
		pter_gcta->fit_reml(grm_file, phen_file, qcovar_file, covar_file, qgxe_file, gxe_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, no_constrain);
	}
	else if(grm_flag || m_grm_flag){
	    if(pca_flag) pter_gcta->pca(grm_file, kp_indi_file, rm_indi_file, grm_cutoff, m_grm_flag, out_pc_num);
	    else if(make_grm_flag) pter_gcta->save_grm(grm_file, kp_indi_file, rm_indi_file, update_sex_file, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, grm_out_bin_flag);
	}
	else throw("Error: no analysis has been launched by the option(s).\n");
	/*      if(genome_flag) pter_gcta->simu_genome(genome_popSize, genome_nSubPOP, genome_nSubSample,
	genome_numPieces, genome_pieceLen, genome_numIndepRegion,
	genome_SNP, genome_rec, genome_mut, genome_mig);
	*/

	delete pter_gcta;
}
