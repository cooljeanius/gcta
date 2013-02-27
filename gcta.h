/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Interface to all the GCTA functions
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#ifndef _GCTA_H
#define _GCTA_H

#ifndef EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#endif

#include "CommFunc.h"
#include "StrFunc.h"
#include "StatFunc.h"
#include <fstream>
#include <iomanip>
#include <bitset>
#include <map>
#include "zfstream.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;
using namespace std;

#ifdef SINGLE_PRECISION
typedef DiagonalMatrix<float, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXf eigenMatrix;
typedef VectorXf eigenVector;
typedef SparseMatrix<float> eigenSparseMat;
typedef DynamicSparseMatrix<float> eigenDynSparseMat;
#else
typedef DiagonalMatrix<double, Dynamic, Dynamic> eigenDiagMat;
typedef MatrixXd eigenMatrix;
typedef VectorXd eigenVector;
typedef SparseMatrix<double> eigenSparseMat;
typedef DynamicSparseMatrix<double> eigenDynSparseMat;
#endif

class gcta
{
public:
    gcta(int autosome_num, double rm_ld_cutoff, string out);
	gcta();
	virtual ~gcta();

    void read_famfile(string famfile);
    void read_bimfile(string bimfile);
    void read_bedfile(string bedfile);
    void read_imp_info_mach(string zinfofile);
    void read_imp_dose_mach(string zdosefile, string kp_indi_file, string rm_indi_file, string blup_indi_file);
    void read_imp_info_beagle(string zinfofile);
    void read_imp_dose_beagle(string zdosefile, string kp_indi_file, string rm_indi_file, string blup_indi_file);
    void update_ref_A(string ref_A_file);
    void update_rsq(string zinfofile);
    void update_freq(string freq);
	void save_freq(bool ssq_flag);
    void extract_snp(string snplistfile);
    void extract_single_snp(string snpname);
    void exclude_snp(string snplistfile);
    void exclude_single_snp(string snpname);
    void extract_chr(int chr_start, int chr_end);
    void filter_snp_maf(double maf);
    void filter_snp_max_maf(double max_maf);
    void filter_dosage_rsq(double rsq_cutoff);
    void keep_indi(string indi_list_file);
    void remove_indi(string indi_list_file);
    void update_sex(string sex_file);
    void read_indi_blup(string blup_indi_file);
	void save_XMat(bool miss_with_mu);

	void paa(string aa_file);
	void ibc(bool ibc_all);
    void make_grm(bool grm_xchr_flag, bool output_bin, int grm_mtd=0);
    void save_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool output_grm_bin);
    void pca(string grm_file, string keep_indi_file, string remove_indi_file, double grm_cutoff, bool merge_grm_flag, int out_pc_num);

    void fit_reml(string grm_file, string phen_file, string qcovar_file, string covar_file, string qGE_file, string GE_file, string keep_indi_file, string remove_indi_file, string sex_file, int mphen, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, int reml_mtd, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, vector<int> drop, bool no_lrt, double prevalence, bool no_constrain, bool mlmassoc=false);
	void blup_snp_geno();
	void blup_snp_dosage();
    void mlm_assoc(string grm_file, string phen_file, string qcovar_file, string covar_file, string keep_indi_file, string remove_indi_file, int mphen, double grm_cutoff, bool m_grm_flag, bool est_fix_eff, int MaxIter, vector<double> reml_priors, double prevalence);

	void read_LD_target_SNPs(string snplistfile);
    void LD_Blocks(int stp, double wind_size, double alpha, bool IncldQ=true, bool save_ram=false);

    void genet_dst(string bfile, string hapmap_genet_map);

    void GWAS_simu(string bfile, int simu_num, string qtl_file, int case_num, int control_num, double hsq, double K, bool output_causal, bool simu_emb_flag);

    void run_massoc_slct(string metafile, int wind_size, double p_cutoff, double collinear, bool joint_only, bool GC, double GC_val, bool actual_geno);
    void run_massoc_cond(string metafile, string snplistfile, int wind_size, double collinear, bool GC, double GC_val, bool actual_geno);
	void run_massoc_sblup(string metafile, int wind_size, double lambda);

    void save_plink();

    void read_IRG_fnames(string snp_info_file, string fname_file, double GC_cutoff);

private:
    void init_keep();
    void init_include();
    void get_rsnp(vector<int> &rsnp);
    void get_rindi(vector<int> &rindi);

    void save_famfile();
    void save_bimfile();
    void save_bedfile();

    void update_bim(vector<int> &rsnp);

    void update_id_map_kp(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
    void update_id_map_rm(const vector<string> &id_list, map<string, int> &id_map, vector<int> &keep);
    void read_snplist(string snplistfile, vector<string> &snplist, string msg="SNPs");
    void read_indi_list(string indi_list_file, vector<string> &indi_list);

    void make_XMat(vector< vector<float> > &X, bool miss_with_mu);
    void make_XMat_SNPs(vector< vector<float> > &X, bool miss_with_mu);
    void std_XMat(vector< vector<float> > &X, vector<double> &sd_SNP, bool grm_xchr_flag, bool divid_by_std=true);
    void makex_eigenVector(int j, eigenVector &x, bool resize=true, bool minus_2p=false);

    //void make_XMat_SNPs(vector< vector<float> > &X, bool miss_with_mu);
    //void std_XMat_SNPs(vector< vector<float> > &X, vector<double> &sd_SNP, bool grm_xchr_flag);
    void calcu_mu(bool ssq_flag=false);
    void mu_func(int j, vector<double> &fac);
    //void ssq_func(int j, vector<double> &fac);
    //void w_func(int j, vector<double> &fac);
	void rm_high_ld();
    void check_autosome();
    void check_chrX();
    void check_sex();

    // grm
    void read_grm_gz(string grm_file, vector<string> &grm_id, bool out_id_log=true, bool read_id_only=false);
    void read_grm_filenames(string merge_grm_file, vector<string> &grm_files, bool out_log=true);
    void merge_grm(string merge_grm_file);
    void rm_cor_indi(double grm_cutoff);
    void adj_grm(double adj_grm_fac);
    void dc(int dosage_compen);
    void manipulate_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag);
    void output_grm_vec(vector< vector<double> > &A, vector< vector<int> > &A_N, bool output_grm_bin);
    void output_grm_MatrixXf(bool output_grm_bin);

    // reml
    void read_phen(string phen_file, vector<string> &phen_ID, vector<string> &phen_buf, int mphen);
    int read_fac(ifstream &ifstrm, vector<string> &ID, vector< vector<string> > &fac);
    int read_covar(string covar_file, vector<string> &covar_ID, vector< vector<string> > &covar, bool qcovar_flag);
    int read_GE(string GE_file, vector<string> &GE_ID, vector< vector<string> > &GE, bool qGE_flag=false);
    bool check_case_control(double &ncase);
    double transform_hsq_L(double P, double K, double hsq);
    int constrain_varcmp(eigenVector &varcmp, double y_Ssq);
    void drop_comp(vector<int> &drop);
    void construct_X(map<string, int> &uni_id_map, bool qcovar_flag, int qcovar_num, vector<string> &qcovar_ID, vector< vector<string> > &qcovar, bool covar_flag, int covar_num, vector<string> &covar_ID, vector< vector<string> > &covar, vector<eigenMatrix> &E_float, eigenMatrix &qE_float);
    void coeff_mat(const vector<string> &vec, eigenMatrix &coeff_mat, string errmsg1, string errmsg2);
    void reml(bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, double prevalence, bool no_constrain, bool no_lrt, bool mlmassoc=false);
    double reml_iteration(double y_Ssq, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp, bool prior_var_flag, bool no_constrain);
    bool calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter);
    bool comput_inverse_logdet_LDLT(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet);
    void bend_A();
    bool bending_eigenval(eigenVector &eval);
    double comput_inverse_logdet_LU(eigenMatrix &Vi, string errmsg);
    double calcu_P(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &P);
	void calcu_Hi(eigenMatrix &P, eigenMatrix &Hi);
	void reml_equation(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp);
	double lgL_reduce_mdl(double y_Ssq, bool no_constrain);
    void em_reml(eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp);
	void ai_reml(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, double dlogL);
	void calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA);
	void calcu_Vp(double &Vp, double &VarVp, int exld_pos, eigenVector &varcmp, eigenMatrix &Hi);
	void calcu_hsq(int i, double Vp, double VarVp, int exld_pos, double &hsq, double &var_hsq, eigenVector &varcmp, eigenMatrix &Hi);
	void output_blup_snp(eigenMatrix &b_SNP);

    // GWAS simulation
    void kosambi();
    int read_QTL_file(string qtl_file, vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, vector<int> &have_eff);
    void output_simu_par(vector<string> &qtl_name, vector<int> &qtl_pos, vector<double> &qtl_eff, double Vp);
    void save_phenfile(vector< vector<double> > &y);
    // not usefully any more
    void GenerCases(string bfile, string qtl_file, int case_num, int control_num, double hsq, double K, bool curr_popu=false, double gnrt=100);

    // estimate LD
	//void Search_LD_Block(int SNP_SmplNum, vector<int> &smpl, int stp, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &v_r2, vector<double> &d, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ);
	void EstLD(vector<int> &smpl, double wind_size, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &max_r2, vector<string> &max_r2_snp, vector<double> &dL, vector<double> &dR, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ);
	void reg(vector<double> &y, vector<double> &x, vector<double> &rst);

    // read raw genotype data (Illumina)
    char flip_allele(char a);
    void read_one_IRG(ofstream &oped, int ind, string IRG_fname, double GC_cutoff);

    // inline functions
    template<typename ElemType>
    void makex(int j, vector<ElemType> &x, bool minus_2p=false)
    {
        int i=0;
        x.resize(_keep.size());
        for(i=0; i<_keep.size(); i++){
            if(!_snp_1[_include[j]][i] || _snp_2[_include[j]][i]){
                if(_allele1[_include[j]]==_ref_A[_include[j]]) x[i]=(_snp_1[_include[j]][i]+_snp_2[_include[j]][i]);
                else x[i]=2.0-(_snp_1[_include[j]][i]+_snp_2[_include[j]][i]);
            }
            else x[i]=_mu[_include[j]];
            if(minus_2p) x[i]-=_mu[_include[j]];
        }
    }

private:
	// read in plink files
	// bim file
	int _autosome_num;
    vector<int> _chr;
    vector<string> _snp_name;
	map<string, int> _snp_name_map;
	vector<double> _genet_dst;
	vector<int> _bp;
	vector<char> _allele1;
	vector<char> _allele2;
	vector<char> _ref_A; // reference allele
	int _snp_num;
	vector<double> _rc_rate;
	vector<int> _include; // initialized in the read_bimfile()

	// fam file
	vector<string> _fid;
	vector<string> _pid;
	map<string, int> _id_map;
	vector<string> _fa_id;
	vector<string> _mo_id;
	vector<int> _sex;
	vector<double> _pheno;
	int _indi_num;
	vector<int> _keep; // initialized in the read_famfile()
	eigenMatrix _varcmp_Py; // BLUP solution to the total genetic effects of individuals

    // bed file
    vector< vector<bool> > _snp_1;
    vector< vector<bool> > _snp_2;

    // dosage file
    bool _dosage_flag;
    vector<double> _dosage_Rsq;
    vector< vector<float> > _geno_dose;

	// QC
	double _rm_ld_cutoff;

    // grm
    eigenMatrix _grm_N;
    eigenMatrix _grm;

    // reml
    int _n;
	int _X_c;
	vector<int> _r_indx;
	vector<int> _r_indx_drop;
	int _reml_max_iter;
	int _reml_mtd;
	int _reml_inv_mtd;
	eigenMatrix _X;
	eigenMatrix _A;
	eigenVector _y;
    eigenMatrix _Vi;
    eigenMatrix _P;
	vector<string> _var_name;
	vector<string> _hsq_name;

    vector<double> _mu;
    string _out;
    bool _save_ram;

    // LD
    vector<string> _ld_target_snp;
};

#endif
