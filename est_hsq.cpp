/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for REML analysis
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_phen(string phen_file, vector<string> &phen_ID, vector<string> &phen_buf, int mphen)
{
    // Read phenotype data
 	ifstream in_phen(phen_file.c_str());
	if(!in_phen) throw("Error: can not open the file ["+phen_file+"] to read.");

    int i=0;
    vector<string> fid, pid, vs_buf;
	string str_buf, fid_buf, pid_buf;
	phen_ID.clear();
	phen_buf.clear();
    cout<<"Reading phenotypes from ["+phen_file+"]."<<endl;
    getline(in_phen, str_buf);
    int phen_num=StrFunc::split_string(str_buf, vs_buf)-2;
    if(phen_num<=0) throw("Error: no phenotype data is found.");
    if(phen_num>1) cout<<"There are "<<phen_num<<" traits specified in the file ["+phen_file+"]."<<endl;
    if(mphen>phen_num){
        stringstream errmsg;
        errmsg<<"Error: can not find the "<<mphen<<"th trait in the file ["+phen_file+"].";
        throw(errmsg.str());
    }
    if(phen_num>1) cout<<"The "<<mphen<<"th trait is included for analysis."<<endl;
    in_phen.seekg(ios::beg);
    mphen--;
	while(in_phen){
	    in_phen>>fid_buf;
		if(in_phen.eof()) break;
		in_phen>>pid_buf;
		for(i=0; i<phen_num; i++){
		    in_phen>>str_buf;
            if(i==mphen && str_buf!="-9" && str_buf!="NA"){
                phen_ID.push_back(fid_buf+":"+pid_buf);
                fid.push_back(fid_buf);
                pid.push_back(pid_buf);
                phen_buf.push_back(str_buf);
            }
		}
	}
	in_phen.close();
    cout<<"Nonmissing phenotypes of "<<phen_buf.size()<<" individuals are included from ["+phen_file+"]."<<endl;

    if(_id_map.empty()){
        _fid=fid;
        _pid=pid;
        _indi_num=_fid.size();
        init_keep();
    }
}

int gcta::read_covar(string covar_file, vector<string> &covar_ID, vector< vector<string> > &covar, bool qcovar_flag)
{
    // Read phenotype data
 	ifstream in_covar(covar_file.c_str());
	if(!in_covar) throw("Error: can not open the file ["+covar_file+"] to read.");

	int i=0, covar_num=0;
	string str_buf, id_buf;
	vector<string> covar_buf, vs_buf;
	covar_ID.clear();
	covar.clear();
    if(qcovar_flag) cout<<"Reading quantative covariates from ["+covar_file+"]."<<endl;
    else cout<<"Reading discrete covariates from ["+covar_file+"]."<<endl;
    covar_num=read_fac(in_covar, covar_ID, covar);
    if(qcovar_flag) cout<<covar_num<<" quantative covarites of "<<covar_ID.size()<<" individuals read from ["+covar_file+"]."<<endl;
    else cout<<covar_num<<" discrete covarites of "<<covar_ID.size()<<" individuals are included from ["+covar_file+"]."<<endl;

    return covar_num;
}

int gcta::read_fac(ifstream &ifstrm, vector<string> &ID, vector< vector<string> > &fac)
{
	int i=0, line=0, fac_num=0, prev_fac_num=0;
	string str_buf, id_buf;
	vector<string> vs_buf;
	while(ifstrm){
	    ifstrm>>str_buf;
		if(ifstrm.eof()) break;
	    id_buf=str_buf;
		ifstrm>>str_buf;
		id_buf+=":"+str_buf;
		getline(ifstrm, str_buf);
		fac_num=StrFunc::split_string(str_buf, vs_buf);
		if(line>0 && fac_num!=prev_fac_num) throw("Error: each row should have the same number of columns.\n"+id_buf+"\t"+str_buf);
		line++;
		prev_fac_num=fac_num;
		bool continue_flag=false;
		for(i=0; i<fac_num; i++){
		    if(vs_buf[i]=="-9" || vs_buf[i]=="NA") continue_flag=true;
		}
		if(continue_flag) continue;
		ID.push_back(id_buf);
		fac.push_back(vs_buf);
	}
	ifstrm.close();
	return fac_num;
}

int gcta::read_GE(string GE_file, vector<string> &GE_ID, vector< vector<string> > &GE, bool qGE_flag)
{
    // Read phenotype data
 	ifstream in_GE(GE_file.c_str());
	if(!in_GE) throw("Error: can not open the file ["+GE_file+"] to read.");

	string str_buf, id_buf;
	vector<string> vs_buf;
	GE_ID.clear();
	GE.clear();
	string env="environmental";
	if(qGE_flag==true) env="continuous "+env;
	else env="categorical "+env;
    cout<<"Reading "<<env<<" factor(s) for the analysis of GE interaction from ["+GE_file+"]."<<endl;
    int GE_num=read_fac(in_GE, GE_ID, GE);
    if(GE_num==0) throw("Error: no "+env+" factor is specified. Please check the format of the file: "+GE_file+".");
    cout<<GE_num<<" "<<env<<" factor(s) for "<<GE_ID.size()<<" individuals are included from ["+GE_file+"]."<<endl;

    return GE_num;
}

void gcta::fit_reml(string grm_file, string phen_file, string qcovar_file, string covar_file, string qGE_file, string GE_file, string keep_indi_file, string remove_indi_file, string sex_file, int mphen, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool m_grm_flag, bool pred_rand_eff, bool est_fix_eff, int reml_mtd, int MaxIter, vector<double> reml_priors, vector<double> reml_priors_var, vector<int> drop, bool no_lrt, double prevalence, bool no_constrain, bool mlmassoc)
{
    _reml_mtd=reml_mtd;
    _reml_max_iter=MaxIter;
    int i=0, j=0, k=0;
    bool grm_flag=(!grm_file.empty());
    bool qcovar_flag=(!qcovar_file.empty());
    bool covar_flag=(!covar_file.empty());
    bool GE_flag=(!GE_file.empty());
    bool qGE_flag=(!qGE_file.empty());
    if(m_grm_flag) grm_flag=false;

    // Read data
    stringstream errmsg;
    int qcovar_num=0, covar_num=0, qE_fac_num=0, E_fac_num=0;
    vector<string> phen_ID, phen_buf, qcovar_ID, covar_ID, qGE_ID, GE_ID, grm_id, grm_files;
    vector< vector<string> > qcovar, covar, GE, qGE; // save individuals by column

    if(grm_flag){
        read_grm_gz(grm_file, grm_id);
        update_id_map_kp(grm_id, _id_map, _keep);
        grm_files.push_back(grm_file);
    }
    else if(m_grm_flag){
        read_grm_filenames(grm_file, grm_files, false);
        for(i=0; i<grm_files.size(); i++){
            read_grm_gz(grm_files[i], grm_id, false, true);
            update_id_map_kp(grm_id, _id_map, _keep);
        }
    }
    read_phen(phen_file, phen_ID, phen_buf, mphen);
    update_id_map_kp(phen_ID, _id_map, _keep);
    if(qcovar_flag){
        qcovar_num=read_covar(qcovar_file, qcovar_ID, qcovar, true);
        update_id_map_kp(qcovar_ID, _id_map, _keep);
    }
    if(covar_flag){
        covar_num=read_covar(covar_file, covar_ID, covar, false);
        update_id_map_kp(covar_ID, _id_map, _keep);
    }
    if(qGE_flag){
        qE_fac_num=read_GE(qGE_file, qGE_ID, qGE, true);
        update_id_map_kp(qGE_ID, _id_map, _keep);
    }
    if(GE_flag){
        E_fac_num=read_GE(GE_file, GE_ID, GE, false);
        update_id_map_kp(GE_ID, _id_map, _keep);
    }
    if(!mlmassoc){
        if(!keep_indi_file.empty()) keep_indi(keep_indi_file);
        if(!remove_indi_file.empty()) remove_indi(remove_indi_file);       
    }
    if(grm_flag){
        if(grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
        if(!sex_file.empty()) update_sex(sex_file);
        if(adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
        if(dosage_compen>-1) dc(dosage_compen);
        _grm_N.resize(1,1);
    }

    vector<string> uni_id;
	map<string, int> uni_id_map;
    map<string, int>::iterator iter;
	for(i=0; i<_keep.size(); i++){
	    uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
	    uni_id_map.insert(pair<string,int>(_fid[_keep[i]]+":"+_pid[_keep[i]], i));
	}
    _n=_keep.size();
    if(_n<1) throw("Error: no individual is in common in the input files.");

    // construct model terms
    _y.setZero(_n);
    for(i=0; i<phen_ID.size(); i++){
        iter=uni_id_map.find(phen_ID[i]);
        if(iter==uni_id_map.end()) continue;
        _y[iter->second]=atof(phen_buf[i].c_str());
    }

    int pos=0;
    _r_indx.clear();
    eigenMatrix A_N(_n, _n);
    vector<int> kp;
    if(grm_flag){
        for(i=0; i<1+qE_fac_num+E_fac_num; i++) _r_indx.push_back(i);
        if(!no_lrt) drop_comp(drop);
        _A=eigenMatrix::Zero(_n, _r_indx.size()*_n);
        if(mlmassoc) StrFunc::match(uni_id, grm_id, kp);
        else kp=_keep;
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++) (_A.block(0,0,_n,_n))(j,i)=(_A.block(0,0,_n,_n))(i,j)=_grm(kp[i],kp[j]);
        }
        pos++;
        _grm.resize(1,1);
    }
    else if(m_grm_flag){
        if(!sex_file.empty()) update_sex(sex_file);
        for(i=0; i<(1+qE_fac_num+E_fac_num)*grm_files.size(); i++) _r_indx.push_back(i);
        if(!no_lrt) drop_comp(drop);
        _A=eigenMatrix::Zero(_n, _r_indx.size()*_n);
        string prev_file=grm_files[0];
        vector<string> prev_grm_id(grm_id);
        cout<<"There are "<<grm_files.size()<<" GRM file names specified in the file ["+grm_file+"]."<<endl;
        for(i=0; i<grm_files.size(); i++, pos++){
            cout<<"Reading the GRM from the "<<i+1<<"th file ..."<<endl;
            read_grm_gz(grm_files[i], grm_id);
            if(adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
            if(dosage_compen>-1) dc(dosage_compen);
            StrFunc::match(uni_id, grm_id, kp);
            int pos_n=pos*_n;
            for(j=0; j<_n; j++){
                for(k=0; k<=j; k++){
                    if(kp[j]>=kp[k]) (_A.block(0,pos_n,_n,_n))(k,j)=(_A.block(0,pos_n,_n,_n))(j,k)=_grm(kp[j],kp[k]);
                    else (_A.block(0,pos_n,_n,_n))(k,j)=(_A.block(0,pos_n,_n,_n))(j,k)=_grm(kp[k],kp[j]);
                }
            }
            prev_file=grm_files[i];
            prev_grm_id=grm_id;
       }
        _grm_N.resize(1,1);
        _grm.resize(1,1);
    }

    // GE interaction
    vector<eigenMatrix> E_float(E_fac_num);
    eigenMatrix qE_float, mbuf;
    if(qGE_flag){
        qE_float.resize(_n, qE_fac_num);
        for(i=0; i<qGE_ID.size(); i++){
            iter=uni_id_map.find(qGE_ID[i]);
            if(iter==uni_id_map.end()) continue;
            for(j=0; j<qE_fac_num; j++) qE_float(iter->second,j)=atof(qGE[i][j].c_str());
        }
        for(j=0; j<qE_fac_num; j++){
            mbuf=((qE_float.block(0,j,_n,1))*(qE_float.block(0,j,_n,1)).transpose());
            for(i=0; i<grm_files.size(); i++, pos++) (_A.block(0,static_cast<Index>(pos)*_n,_n,_n))=(_A.block(0,static_cast<Index>(i)*_n,_n,_n)).array()*mbuf.array();
        }
    }
    if(GE_flag){
        vector< vector<string> > E_str(E_fac_num);
        for(i=0; i<E_fac_num; i++) E_str[i].resize(_n);
        for(i=0; i<GE_ID.size(); i++){
            iter=uni_id_map.find(GE_ID[i]);
            if(iter!=uni_id_map.end()){
                for(j=0; j<E_fac_num; j++) E_str[j][iter->second]=GE[i][j];
            }
        }
        for(j=0; j<E_fac_num; j++){
            stringstream errmsg;
            errmsg<<"Error: too many classes for the "<<j+1<<"th environmental factor. \nPlease make sure you input a discrete variable as the environmental factor.";
            string errmsg1=errmsg.str();
            errmsg.str("");
            errmsg<<"Error: the "<<j+1<<"th envronmental factor has only one class.";
            string errmsg2=errmsg.str();
            coeff_mat(E_str[j], E_float[j], errmsg1, errmsg2);
            mbuf=((E_float[j])*(E_float[j]).transpose());
            for(i=0; i<grm_files.size(); i++, pos++) (_A.block(0,static_cast<Index>(pos)*_n,_n,_n))=(_A.block(0,static_cast<Index>(i)*_n,_n,_n)).array()*mbuf.array();
        }
    }

    // construct X matrix
    construct_X(uni_id_map, qcovar_flag, qcovar_num, qcovar_ID, qcovar, covar_flag, covar_num, covar_ID, covar, E_float, qE_float);

    // names of variance component
    for(i=0; i<grm_files.size(); i++){
        stringstream strstrm;
        if(grm_files.size()==1) strstrm<<"";
        else strstrm<<i+1;
        _var_name.push_back("V(G"+strstrm.str()+")");
        _hsq_name.push_back("V(G"+strstrm.str()+")/Vp");
    }
    for(j=0; j<qE_fac_num; j++){
        for(i=0; i<grm_files.size(); i++){
            stringstream strstrm1,strstrm2;
            if(grm_files.size()==1) strstrm1<<"";
            else  strstrm1<<i+1;
            if(qE_fac_num==1) strstrm2<<"";
            else strstrm2<<j+1;
            _var_name.push_back("V(G"+strstrm1.str()+"xqE"+strstrm2.str()+")");
            _hsq_name.push_back("V(G"+strstrm1.str()+"xqE"+strstrm2.str()+")"+"/Vp");
        }
    }
    for(j=0; j<E_fac_num; j++){
        for(i=0; i<grm_files.size(); i++){
            stringstream strstrm1,strstrm2;
            if(grm_files.size()==1) strstrm1<<"";
            else  strstrm1<<i+1;
            if(E_fac_num==1) strstrm2<<"";
            else strstrm2<<j+1;
            _var_name.push_back("V(G"+strstrm1.str()+"xE"+strstrm2.str()+")");
            _hsq_name.push_back("V(G"+strstrm1.str()+"xE"+strstrm2.str()+")"+"/Vp");
        }
    }
    _var_name.push_back("V(e)");

    // run REML algorithm
    cout<<_n<<" individuals are in common in these files."<<endl;
	reml(pred_rand_eff, est_fix_eff, reml_priors, reml_priors_var, prevalence, no_constrain, no_lrt, mlmassoc);
}

void gcta::drop_comp(vector<int> &drop)
{
    int i=0;
    stringstream errmsg;
    _r_indx_drop=_r_indx;
    stable_sort(drop.begin(), drop.end());
    drop.erase(unique(drop.begin(), drop.end()), drop.end());
    for(i=drop.size()-1; i>=0; i--){
        if(drop[i]<1 || drop[i]>_r_indx.size()){
            errmsg<<"Error: there are "<<_r_indx.size()<<" genetic variance component in the model and can't drop the "<<drop[i]<<"-th component.";
            throw(errmsg.str());
        }
        _r_indx_drop.erase(_r_indx_drop.begin()+drop[i]-1);
    }
    if(_r_indx.size()==_r_indx_drop.size()) throw("Erorr: no component has been dropped from the model. Please check the --reml-lrt option.");
}

void gcta::construct_X(map<string, int> &uni_id_map, bool qcovar_flag, int qcovar_num, vector<string> &qcovar_ID, vector< vector<string> > &qcovar, bool covar_flag, int covar_num, vector<string> &covar_ID, vector< vector<string> > &covar, vector<eigenMatrix> &E_float, eigenMatrix &qE_float)
{
    int i=0, j=0;
    map<string, int>::iterator iter;
    stringstream errmsg;

    _X_c=1;
    // quantitative covariates
    eigenMatrix X_q;
    if(qcovar_flag){
        X_q.resize(_n, qcovar_num);
        for(i=0; i<qcovar_ID.size(); i++){
            iter=uni_id_map.find(qcovar_ID[i]);
            if(iter==uni_id_map.end()) continue;
            for(j=0; j<qcovar_num; j++) X_q(iter->second,j)=atof(qcovar[i][j].c_str());
        }
        if(qcovar_num==0) throw("Error: no quantitative covariate is found.");
        cout<<qcovar_num<<" quantitative variable(s) included as covariate(s)."<<endl;
        _X_c+=qcovar_num;
    }
    // discrete covariates
    vector<eigenMatrix> X_d;
    if(covar_flag){
        vector< vector<string> > covar_tmp(covar_num);
        for(i=0; i<covar_num; i++) covar_tmp[i].resize(_n);
        for(i=0; i<covar_ID.size(); i++){
            iter=uni_id_map.find(covar_ID[i]);
            if(iter==uni_id_map.end()) continue;
            for(j=0; j<covar_num; j++) covar_tmp[j][iter->second]=covar[i][j];
        }
        cout<<covar_num<<" discrete variable(s) included as covariate(s)."<<endl;
        X_d.resize(covar_num);
        for(i=0; i<covar_num; i++){
            stringstream errmsg;
            errmsg<<"Error: too many classes for the "<<i+1<<"th discrete variable. \nPlease use the --qcovariate if it is a quantitative covariate.";
            string errmsg1=errmsg.str();
            errmsg.str("");
            errmsg<<"Error: the "<<i+1<<"th discrete variable has only one class.";
            string errmsg2=errmsg.str();
            coeff_mat(covar_tmp[i], X_d[i], errmsg1, errmsg2);
            _X_c+=(X_d[i]).cols()-1;
        }
    }

    // E factor
    _X_c+=qE_float.cols();
    for(i=0; i<E_float.size(); i++) _X_c+=(E_float[i]).cols()-1;

    // Construct _X
    int col=0;
    _X.resize(_n, _X_c);
    _X.block(0, col, _n, 1)=eigenMatrix::Ones(_n, 1);
    col++;
    if(qcovar_flag){
        _X.block(0, col, _n, X_q.cols())=X_q;
        col+=X_q.cols();
    }
    for(i=0; i<X_d.size(); i++){
        _X.block(0, col, _n, (X_d[i]).cols()-1)=(X_d[i]).block(0,1,_n,(X_d[i]).cols()-1);
        col+=(X_d[i]).cols()-1;
    }
    if(qE_float.cols()>0){
        _X.block(0, col, _n, qE_float.cols())=qE_float;
        col+=qE_float.cols();
    }
    for(i=0; i<E_float.size(); i++){
        _X.block(0, col, _n, (E_float[i]).cols()-1)=(E_float[i]).block(0,1,_n,(E_float[i]).cols()-1);
        col+=(E_float[i]).cols()-1;
    }
}

void gcta::coeff_mat(const vector<string> &vec, eigenMatrix &coeff_mat, string errmsg1, string errmsg2)
{
	vector<string> value(vec);
	stable_sort(value.begin(), value.end());
	value.erase(unique(value.begin(), value.end()), value.end());
	if(value.size()>0.5*vec.size()) throw(errmsg1); // throw("Error: too many classes for the envronmental factor. \nPlease make sure you input a discrete variable as the environmental factor.");
	if(value.size()==1) throw(errmsg2); //throw("Error: the envronmental factor should has more than one classes.");

	int i=0, j=0, row_num=vec.size(), column_num=value.size();
	map<string, int> val_map;
	for(i=0; i<value.size(); i++) val_map.insert(pair<string, int>(value[i], i));

    coeff_mat.resize(row_num, column_num);
	coeff_mat.setZero(row_num, column_num);
	map<string, int>::iterator iter;
	for(i=0; i<row_num; i++){
		iter=val_map.find(vec[i]);
		coeff_mat(i,iter->second)=1.0;
	}
}

bool gcta::check_case_control(double &ncase)
{
    double case_num=0.0;
	vector<double> value(_n);
	for(int i=0; i<_n; i++) value[i]=_y(i);
	stable_sort(value.begin(), value.end());
	value.erase(unique(value.begin(), value.end()), value.end());
	if(value.size()==2){
	    if(CommFunc::FloatEqual(value[0],0.0) && CommFunc::FloatEqual(value[1],1.0)) case_num=_y.sum();
	    else if(CommFunc::FloatEqual(value[0],1.0) && CommFunc::FloatEqual(value[1],2.0)) case_num=(_y.sum()-_n);
        cout<<"Assuming a disease phenotype for a case-control study ";
        cout<<"("<<(int)case_num<<" cases and "<<(int)(_n-case_num)<<" controls)."<<endl;
        ncase=case_num/_n;
        return true;
    }
	else if(value.size()<2) throw("Error: invalid phenotype. Please check the phenotype file.");
	return false;
}

double gcta::transform_hsq_L(double P, double K, double hsq)
{
    double t=StatFunc::qnorm(1.0-K);
    double z=StatFunc::dnorm(t);
    double C=(K*(1-K)/(z*z))*(K*(1-K)/(P*(1-P)));
    return(hsq*C);
}

int gcta::constrain_varcmp(eigenVector &varcmp, double y_Ssq)
{
    double delta=0.0;
    int i=0, num=0;
    vector<int> constrain(_r_indx.size()+1);
    for(i=0; i<_r_indx.size()+1; i++){
        if(varcmp[i]<0){
            delta+=y_Ssq*1e-6-varcmp[i];
            varcmp[i]=y_Ssq*1e-6;
            constrain[i]=1;
            num++;
        }
    }
    delta/=(_r_indx.size()+1-num);
    for(i=0; i<_r_indx.size()+1; i++){
        if(constrain[i]<1 && varcmp[i]>delta) varcmp[i]-=delta;
    }

    return num;
}

void gcta::reml(bool pred_rand_eff, bool est_fix_eff, vector<double> &reml_priors, vector<double> &reml_priors_var, double prevalence, bool no_constrain, bool no_lrt, bool mlmassoc)
{
	int i=0, j=0, k=0;

    // case-control
    double ncase=0.0;
    bool flag_cc=check_case_control(ncase);
    if(flag_cc){
        if(mlmassoc) throw("Error: the option --mlm-assoc is valid for the quantitative trait only.");
        if(prevalence<-1) cout<<"Warning: please specify the disease prevalence by the option --prevalence so that GCTA can transform the variance explained to the underlying liability scale."<<endl;
    }

	// Initialize variance component
	// 0: AI; 1: REML equation; 2: EM
    stringstream errmsg;
	double d_buf=0.0, y_mean=_y.mean();
	eigenVector y(_y);
	for(i=0; i<_n; i++) y(i)-=y_mean;
	double y_Ssq=y.squaredNorm()/(_n-1.0);
	if(!(fabs(y_Ssq)<1e30)) throw("Error: the phenotypic variance is infinite. Please check the missing data in your phenotype file. Missing values should be represented by \"NA\" or \"-9\".");
	bool reml_priors_flag=!reml_priors.empty(), reml_priors_var_flag=!reml_priors_var.empty();
	if(reml_priors_flag && reml_priors.size()<_r_indx.size()){
	    errmsg<<"Error: in option --reml-priors. There are "<<_r_indx.size()+1<<" variance components. At least "<<_r_indx.size()<<" prior values should be specified.";
	    throw(errmsg.str());
	}
	if(reml_priors_var_flag && reml_priors_var.size()<_r_indx.size()){
	    errmsg<<"Error: in option --reml-priors-var. There are "<<_r_indx.size()+1<<" variance components. At least "<<_r_indx.size()<<" prior values should be specified.";
	    throw(errmsg.str());
	}

	cout<<"\nPerforming REML analysis ... (NOTE: may take hours depending on sample size)."<<endl;
	if(_n<10) throw("Error: sample size is too small.");
	cout<<_n<<" observations, "<<_X_c<<" fixed effect(s), and "<<_r_indx.size()+1<<" variance component(s)(including residual variance)."<<endl;
	eigenMatrix Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(_r_indx.size()+1, _r_indx.size()+1);
    eigenVector Py(_n), varcmp(_r_indx.size()+1);
	if(reml_priors_var_flag){
	    for(i=0; i<_r_indx.size(); i++) varcmp[i]=reml_priors_var[i];
	    if(varcmp[_r_indx.size()]<1e-30) varcmp[i]=y_Ssq-varcmp.sum();
	}
	else if(reml_priors_flag){
	    for(i=0; i<_r_indx.size(); i++) { varcmp[i]=reml_priors[i]*y_Ssq; d_buf+=reml_priors[i]; }
	    varcmp[_r_indx.size()]=(1.0-d_buf)*y_Ssq;
	}
	else varcmp.setConstant(y_Ssq/(_r_indx.size()+1));
	double lgL=reml_iteration(y_Ssq, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, reml_priors_var_flag|reml_priors_flag, no_constrain);
    eigenMatrix u;
    if(pred_rand_eff){
        u.resize(_n, _r_indx.size());
        for(i=0; i<_r_indx.size(); i++) (u.col(i))=(((_A.block(0, static_cast<Index>(_r_indx[i]) * _n, _n, _n))*Py)*varcmp[i]);
    }
	eigenVector b;
	if(est_fix_eff){
	    b.resize(_X_c);
	    b=Xt_Vi_X_i*(Vi_X.transpose()*_y);
	}
    // calculate Hsq and SE
	double Vp=0.0, VarVp=0.0, Vp_f=0.0, VarVp_f=0.0;
	vector<double> Hsq(_r_indx.size()), VarHsq(_r_indx.size());
	calcu_Vp(Vp, VarVp, -1, varcmp, Hi);
	for(i=0; i<_r_indx.size(); i++) calcu_hsq(i, Vp, VarVp, -1, Hsq[i], VarHsq[i], varcmp, Hi);

	// calculate the logL for a reduce model
	double lgL_rdu_mdl=0.0, LRT=0.0;
	if(!no_lrt){
	    lgL_rdu_mdl=lgL_reduce_mdl(y_Ssq, no_constrain);
	    LRT=2.0*(lgL-lgL_rdu_mdl);
	    if(LRT<0.0) LRT=0.0;
	}

	// output results
	cout<<"\nSummary result of REML analysis:"<<endl;
	cout<<"Source\tVariance\tSE"<<setiosflags(ios::fixed)<<setprecision(6)<<endl;
	for(i=0; i<_r_indx.size()+1; i++) cout<<_var_name[i]<<"\t"<<varcmp[i]<<"\t"<<sqrt(Hi(i,i))<<endl;
	cout<<"Vp\t"<<Vp<<"\t"<<sqrt(VarVp)<<endl;
	for(i=0; i<_r_indx.size(); i++) cout<<_hsq_name[i]<<"\t"<<Hsq[i]<<"\t"<<sqrt(VarHsq[i])<<endl;
	if(flag_cc && prevalence>-1){
        cout<<"The estimate of variance explained on the observed scale is transformed to that on the underlying scale:"<<endl;
        cout<<"(Proportion of cases in the sample = "<<ncase<<"; User-specified disease prevalence = "<<prevalence<<")"<<endl;
	    for(i=0; i<_r_indx.size(); i++) cout<<_hsq_name[i]<<"_L\t"<<transform_hsq_L(ncase, prevalence, Hsq[i])<<"\t"<<transform_hsq_L(ncase, prevalence, sqrt(VarHsq[i]))<<endl;
	}
    if(mlmassoc) return;
	cout<<"\nCovariance/Variance/Correlation Matrix:"<<endl;
	for(i=0; i<_r_indx.size()+1; i++){
	    for(j=0; j<=i; j++) cout<<setiosflags(ios::scientific)<<Hi(i,j)<<"\t";
	    cout<<endl;
	}
	if(est_fix_eff){
        cout<<"Estimate"<<(_X_c>1?"s":"")<<"of fixed effect"<<(_X_c>1?"s":"")<<":"<<endl;
        cout<<"\nSource\tEstimate\tSE"<<endl;
        for(i=0; i<_X_c; i++){
            if(i==0) cout<<"mean\t";
            else cout<<"X_"<<i+1<<"\t";
            cout<<setiosflags(ios::fixed)<<b[i]<<"\t"<<sqrt(Xt_Vi_X_i(i,i))<<endl;          
        }
	}

	// save summary result into a file
	string reml_rst_file=_out+".hsq";
	ofstream o_reml(reml_rst_file.c_str());
	o_reml<<"Source\tVariance\tSE"<<setiosflags(ios::fixed)<<setprecision(6)<<endl;
	for(i=0; i<_r_indx.size()+1; i++) o_reml<<_var_name[i]<<"\t"<<varcmp[i]<<"\t"<<sqrt(Hi(i,i))<<endl;
	o_reml<<"Vp\t"<<Vp<<"\t"<<sqrt(VarVp)<<endl;
	for(i=0; i<_r_indx.size(); i++) o_reml<<_hsq_name[i]<<"\t"<<Hsq[i]<<"\t"<<sqrt(VarHsq[i])<<endl;
	if(flag_cc && prevalence>-1){
        o_reml<<"The estimate of variance explained on the observed scale is transformed to that on the underlying scale:"<<endl;
        o_reml<<"(Proportion of cases in the sample = "<<ncase<<"; User-specified disease prevalence = "<<prevalence<<")"<<endl;
	    for(i=0; i<_r_indx.size(); i++) o_reml<<_hsq_name[i]<<"_L\t"<<transform_hsq_L(ncase, prevalence, Hsq[i])<<"\t"<<transform_hsq_L(ncase, prevalence, sqrt(VarHsq[i]))<<endl;
	}
	o_reml<<"logL\t"<<setprecision(3)<<lgL<<endl;
	if(!no_lrt && _r_indx.size()>0){
        o_reml<<"logL0\t"<<setprecision(3)<<lgL_rdu_mdl<<endl;
        o_reml<<"LRT\t"<<setprecision(3)<<LRT<<endl;
        o_reml<<"Pval\t"<<setiosflags(ios::scientific)<<0.5*StatFunc::chi_prob(_r_indx.size()-_r_indx_drop.size(), LRT)<<setiosflags(ios::fixed)<<endl;
	}
	o_reml<<"n\t"<<_n<<endl;
	if(est_fix_eff){
        o_reml<<"\nFix_eff\tSE"<<endl;
        for(i=0; i<_X_c; i++) o_reml<<setprecision(6)<<b[i]<<"\t"<<sqrt(Xt_Vi_X_i(i,i))<<endl;
        o_reml.close();
	}
	cout<<"\nSummary result of REML analysis has been saved in the file ["+reml_rst_file+"]."<<endl;

	// save random effect to a file
	if(pred_rand_eff){
        string rand_eff_file=_out+".indi.blp";
        ofstream o_rand_eff(rand_eff_file.c_str());
        for(i=0; i<_keep.size(); i++){
            o_rand_eff<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t";
            for(j=0; j<_r_indx.size(); j++) o_rand_eff<<setprecision(6)<<Py[i]*varcmp[j]<<"\t"<<u(i,j)<<"\t";
            o_rand_eff<<endl;
        }
        cout<<"\nBLUP of the genetic effects for "<<_keep.size()<<" individuals has been saved in the file ["+rand_eff_file+"]."<<endl;
	}
}

double gcta::lgL_reduce_mdl(double y_Ssq, bool no_constrain)
{
    if(_r_indx.size()==0) return 0;
    bool multi_comp=(_r_indx.size()-_r_indx_drop.size()>1);
    cout<<"\nCalculating the logLikelihood for the reduced model ...\n(variance component"<<(multi_comp?"s ":" ");
    for(int i=0; i<_r_indx.size(); i++){
        if(find(_r_indx_drop.begin(), _r_indx_drop.end(), _r_indx[i])==_r_indx_drop.end()) cout<<_r_indx[i]+1<<" ";
    }
    cout<<(multi_comp?"are":"is")<<" dropped from the model)"<<endl;
    vector<int> vi_buf(_r_indx);
    _r_indx=_r_indx_drop;
	eigenMatrix Vi_X(_n, _X_c), Xt_Vi_X_i(_X_c, _X_c), Hi(_r_indx.size()+1, _r_indx.size()+1);
    eigenVector Py(_n);
    eigenVector varcmp(_r_indx.size()+1);
    varcmp.setConstant(y_Ssq/(_r_indx.size()+1));
    double lgL=reml_iteration(y_Ssq, Vi_X, Xt_Vi_X_i, Hi, Py, varcmp, false, no_constrain);
    _r_indx=vi_buf;
	return(lgL);
}

double gcta::reml_iteration(double y_Ssq, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp, bool prior_var_flag, bool no_constrain)
{
    const char *mtd_str[3]={"AI-REML algorithm", "REML equation ...", "EM-REML algorithm ..."};
    int i=0, constrain_num=0, iter=0, reml_mtd_tmp=_reml_mtd;
    double logdet=0.0, logdet_Xt_Vi_X=0.0, prev_lgL=-1e20, lgL=-1e20, dlogL=1000.0;
    eigenVector prev_varcmp(varcmp), varcomp_init(varcmp);
	_Vi.resize(_n,_n);
    _P(_n,_n);
	for(iter=0; iter<_reml_max_iter; iter++){
	    if(iter==0){
	        prev_varcmp=varcomp_init;
	        if(prior_var_flag){
	            cout<<"User-specified prior values of variance components: "<<varcmp.transpose()<<endl;
	            //continue;
	        }
	        else{
	            _reml_mtd=2;
	            cout<<"Calculating prior values of variance components by EM-REML ..."<<endl;
	        }
	    }
	    if(iter==1){
            _reml_mtd=reml_mtd_tmp;
	        cout<<"Running "<<mtd_str[_reml_mtd]<<" ..."<<"\nIter.\tlogL\t";
            for(i=0; i<_r_indx.size(); i++) cout<<_var_name[_r_indx[i]]<<"\t";
            cout<<_var_name[_var_name.size()-1]<<endl;
        }
		if(!calcu_Vi(_Vi, prev_varcmp, logdet, iter)) continue; // Calculate Vi
		logdet_Xt_Vi_X=calcu_P(_Vi, Vi_X, Xt_Vi_X_i, _P); // Calculate P
		if(_reml_mtd==0) ai_reml(_P, Hi, Py, prev_varcmp, varcmp, dlogL);
		else if(_reml_mtd==1) reml_equation(_P, Hi, Py, varcmp);
		else if(_reml_mtd==2) em_reml(_P, Py, prev_varcmp, varcmp);
		lgL=-0.5*(logdet_Xt_Vi_X+logdet+(_y.transpose()*Py)(0,0));

        // output log
        if(!no_constrain) constrain_num=constrain_varcmp(varcmp, y_Ssq);
        if(iter>0){
            cout<<iter<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<lgL<<"\t";
            for(i=0; i<_r_indx.size()+1; i++) cout<<setprecision(5)<<varcmp[i]<<"\t";
            if(constrain_num>0) cout<<"("<<constrain_num<<" component(s) constrained)"<<endl;
            else cout<<endl;
        }
        else{
            if(!prior_var_flag) cout<<"Prior values updated from EM-REML: "<<varcmp.transpose()<<endl; 
            cout<<"logL: "<<lgL<<endl;
        }
        if(constrain_num*2>_r_indx.size()+1) throw("Error: analysis stopped because more than half of the variance components are constrained. The result would be unreliable.\n Please have a try to add the option --reml-no-constrain.");

		// convergence
		dlogL=lgL-prev_lgL;
		if((varcmp-prev_varcmp).squaredNorm()/varcmp.squaredNorm()<1e-8 && (fabs(dlogL)<1e-4 || (fabs(dlogL)<1e-2 && dlogL<0))){
			if(_reml_mtd==2){ calcu_Hi(_P, Hi); Hi=2*Hi; } // for calculation of SE
			break;
		}
		prev_varcmp=varcmp; prev_lgL=lgL;
	}
	if(iter==_reml_max_iter){
        stringstream errmsg;
        errmsg<<"Error: Log-likelihood not converged (stop after "<<_reml_max_iter<<" iteractions). \nYou can specify the option --reml-maxit to allow for more iterations."<<endl;
        throw(errmsg.str());
    }
	else cout<<"Log-likelihood ratio converged."<<endl;

	return lgL;
}

void gcta::calcu_Vp(double &Vp, double &VarVp, int exld_pos, eigenVector &varcmp, eigenMatrix &Hi)
{
	Vp=0.0; VarVp=0.0;
	for(int i=0; i<_r_indx.size()+1;  i++){
	    if(i==exld_pos) continue;
	    Vp+=varcmp[i];
	    for(int j=0; j<_r_indx.size()+1; j++){
	        if(j==exld_pos) continue;
	        VarVp+=Hi(i,j);
	    }
	}
}

void gcta::calcu_hsq(int i, double Vp, double VarVp, int exld_pos, double &hsq, double &var_hsq, eigenVector &varcmp, eigenMatrix &Hi)
{
	double V1=varcmp[i], VarV1=Hi(i,i), Cov12=0.0;
	for(int j=0; j<_r_indx.size()+1; j++){
	    if(j==exld_pos) continue;
	    Cov12+=Hi(i,j);
	}
	hsq=V1/Vp;
	var_hsq=(V1/Vp)*(V1/Vp)*(VarV1/(V1*V1)+VarVp/(Vp*Vp)-(2*Cov12)/(V1*Vp));
}

bool gcta::calcu_Vi(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet, int &iter)
{
    int i=0, j=0, k=0;
    string errmsg="\nError: the V (variance-covariance) matrix is not invertible.";

    Vi.setZero();
    if(_r_indx.size()==0){
        Vi.diagonal()=eigenVector::Constant(_n, 1.0/prev_varcmp[_r_indx.size()]);
        logdet=_n*log(prev_varcmp[_r_indx.size()]);
    }
    else{
        Vi.diagonal()=eigenVector::Constant(_n, prev_varcmp[_r_indx.size()]);
        for(i=0; i<_r_indx.size(); i++) Vi+=(_A.block(0, static_cast<Index>(_r_indx[i]) * _n, _n, _n)) * prev_varcmp[i];
        if(!comput_inverse_logdet_LDLT(Vi, prev_varcmp, logdet)){
            cout<<"Warning: the variance-covaraince matrix V is negative-definite."<<endl;
            bend_A();
            iter=-1;
            cout<<"Restarting iterations ..."<<endl;
            return false;
        }
    }
    return true;
}

void gcta::bend_A()
{
    cout<<"Bending the GRM(s) to be positive-definite (may take a long time if there are multiple GRMs)..."<<endl;
    int i=0;
    for(i=0; i<_r_indx.size(); i++){
        SelfAdjointEigenSolver<eigenMatrix> eigensolver(_A.block(0,_r_indx[i]*_n,_n,_n));
        eigenVector eval=eigensolver.eigenvalues();
        if(bending_eigenval(eval)){
            _A.block(0,_r_indx[i]*_n,_n,_n)=eigensolver.eigenvectors()*eigenDiagMat(eval)*eigensolver.eigenvectors().transpose();
            cout<<"Bending the "<<i+1<<"th GRM completed."<<endl;
        }
    }
}

bool gcta::bending_eigenval(eigenVector &eval)
{
    int j=0;
    double eval_m=eval.mean();
    if(eval.minCoeff()>0.0) return false;
    double S=0.0, P=0.0;
    for(j=0; j<eval.size(); j++){
        if(eval[j]>=0) continue;
        S+=eval[j];
        P=-eval[j];
    }
    double W=S*S*100.0+1;
    for(j=0; j<eval.size(); j++){
        if(eval[j]>=0) continue;
        eval[j]=P*(S-eval[j])*(S-eval[j])/W;
    }
    eval*=eval_m/eval.mean();
    return true;
}

bool gcta::comput_inverse_logdet_LDLT(eigenMatrix &Vi, eigenVector &prev_varcmp, double &logdet)
{
    int i=0, n=Vi.cols();
    LDLT<eigenMatrix> ldlt(Vi);
    eigenVector d=ldlt.vectorD();
    if(d.minCoeff()<0){
        if(prev_varcmp.minCoeff()>0) return false;
        else throw("Error: the matrix V becomes negative-definite because of one of the variance component is negative.\nPlease re-run the analysis without the --reml-no-constrain option.");
    }
    logdet=0.0;
    for(i=0; i<n; i++) logdet+=log(d[i]);
    Vi.setIdentity();
    ldlt.solveInPlace(Vi);
    return true;
}

double gcta::comput_inverse_logdet_LU(eigenMatrix &Vi, string errmsg)
{
    double logdet=0.0;
    int n=Vi.cols();

    FullPivLU<eigenMatrix> lu(Vi);
    eigenVector u=lu.matrixLU().diagonal();
    if(!lu.isInvertible()) throw(errmsg);
    for(int i=0; i<n; i++) logdet+=log(fabs(u[i]));
    Vi=lu.inverse();
    return logdet;
}

double gcta::calcu_P(eigenMatrix &Vi, eigenMatrix &Vi_X, eigenMatrix &Xt_Vi_X_i, eigenMatrix &P)
{
    Vi_X=Vi*_X;
	Xt_Vi_X_i=_X.transpose()*Vi_X;
	double logdet_Xt_Vi_X=comput_inverse_logdet_LU(Xt_Vi_X_i, "\nError: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).");
	P=Vi-Vi_X*Xt_Vi_X_i*Vi_X.transpose();
	return logdet_Xt_Vi_X;
}

// input P, calculate PA and Hi
void gcta::calcu_Hi(eigenMatrix &P, eigenMatrix &Hi)
{
    int i=0, j=0, k=0, l=0;
    double d_buf=0.0;

	// Calculate PA
    vector<eigenMatrix> PA(_r_indx.size());
	for(i=0; i<_r_indx.size(); i++){
	    (PA[i]).resize(_n, _n);
	    (PA[i])=P*(_A.block(0,_r_indx[i]*_n,_n,_n));
	}

	// Calculate Hi
	for(i=0; i<_r_indx.size(); i++){
		for(j=0; j<=i; j++){
			d_buf=0.0;
			for(k=0; k<_n; k++){
				for(l=0; l<_n; l++) d_buf+=(PA[i])(k,l)*(PA[j])(l,k);
			}
			Hi(i,j)=Hi(j,i)=d_buf;
		}
	}
	for(j=0; j<_r_indx.size(); j++){
		d_buf=0.0;
		for(k=0; k<_n; k++){
			for(l=0; l<_n; l++) d_buf+=P(k,l)*(PA[j])(l,k);
		}
		Hi(j,_r_indx.size())=Hi(_r_indx.size(),j)=d_buf;
	}
	d_buf=0.0;
	for(k=0; k<_n; k++){
		for(l=0; l<k; l++) d_buf+=P(k,l)*P(k,l);
		for(l=k; l<_n; l++) d_buf+=P(l,k)*P(l,k);
	}
	Hi(_r_indx.size(),_r_indx.size())=d_buf;
	comput_inverse_logdet_LU(Hi, "Error: information matrix is not invertible.");
}

// use REML equation to estimate variance component
// input P, calculate PA, H, R and varcmp
void gcta::reml_equation(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &varcmp)
{
	// Calculate Hi
	calcu_Hi(P, Hi);

	// Calculate R
	Py=P*_y;
	eigenVector R(_r_indx.size()+1);
	for(int i=0; i<_r_indx.size(); i++) R(i)=(Py.transpose()*(_A.block(0,_r_indx[i]*_n,_n,_n))*Py)(0,0);
	R(_r_indx.size())=Py.squaredNorm();

	// Calculate variance component
	varcmp=Hi*R;
	Hi=2*Hi; // for calculation of SE
}

// use REML equation to estimate variance component
// input P, calculate PA, H, R and varcmp
void gcta::ai_reml(eigenMatrix &P, eigenMatrix &Hi, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp, double dlogL)
{
    int i=0, j=0;

 	Py=P*_y;
	eigenVector cvec(_n);
	eigenMatrix APy(_n, _r_indx.size());
	for(i=0; i<_r_indx.size(); i++) (APy.col(i))=(_A.block(0,_r_indx[i]*_n,_n,_n))*Py;

	// Calculate Hi
	eigenVector R(_r_indx.size()+1);
	for(i=0; i<_r_indx.size(); i++){
		R(i)=(Py.transpose()*(APy.col(i)))(0,0);
		cvec=P*(APy.col(i));
		Hi(i,i)=((APy.col(i)).transpose()*cvec)(0,0);
		for(j=0; j<i; j++) Hi(j,i)=Hi(i,j)=((APy.col(j)).transpose()*cvec)(0,0);
	}
	cvec=P*Py;
	for(j=0; j<_r_indx.size(); j++) Hi(j,_r_indx.size())=Hi(_r_indx.size(),j)=((APy.col(j)).transpose()*cvec)(0,0);
	R(_r_indx.size())=(Py.transpose()*Py)(0,0);
	Hi(_r_indx.size(),_r_indx.size())=(Py.transpose()*cvec)(0,0);
	Hi=0.5*Hi;

	// Calcualte tr(PA) and dL
    eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA);
	R=-0.5*(tr_PA-R);

	// Calculate variance component
	comput_inverse_logdet_LU(Hi, "Error: AI matrix is not invertible.");
	eigenVector delta(_r_indx.size()+1);
	delta=Hi*R;
	if(dlogL>1.0) varcmp=prev_varcmp+0.316*delta;
	else varcmp=prev_varcmp+delta;
}

// input P, calculate varcmp
void gcta::em_reml(eigenMatrix &P, eigenVector &Py, eigenVector &prev_varcmp, eigenVector &varcmp)
{
	int i=0;

 	// Calculate trace(PA)
 	eigenVector tr_PA;
	calcu_tr_PA(P, tr_PA);

	// Calculate R
	Py=P*_y;
	eigenVector R(_r_indx.size()+1);
	for(i=0; i<_r_indx.size(); i++) R(i)=(Py.transpose()*(_A.block(0,_r_indx[i]*_n,_n,_n))*Py)(0,0);
	R(_r_indx.size())=Py.squaredNorm();

	// Calculate variance component
	for(i=0; i<_r_indx.size(); i++) varcmp(i)=(prev_varcmp(i)*_n-prev_varcmp(i)*prev_varcmp(i)*tr_PA(i)+prev_varcmp(i)*prev_varcmp(i)*R(i))/_n;
	int j=_r_indx.size();
	varcmp(j)=(prev_varcmp(j)*_n-prev_varcmp(j)*prev_varcmp(j)*tr_PA(j)+prev_varcmp(j)*prev_varcmp(j)*R(j))/_n;
}

// input P, calculate tr(PA)
void gcta::calcu_tr_PA(eigenMatrix &P, eigenVector &tr_PA)
{
	int i=0, k=0, l=0;
    double d_buf=0.0;

	// Calculate trace(PA)
	tr_PA.resize(_r_indx.size()+1);
	for(i=0; i<_r_indx.size(); i++){
		d_buf=0.0;
		for(k=0; k<_n; k++){
			for(l=0; l<_n; l++) d_buf+=P(k,l)*(_A.block(0,_r_indx[i]*_n,_n,_n))(k,l);
		}
		tr_PA(i)=d_buf;
	}
	tr_PA(_r_indx.size())=P.trace();
}

// blue estimate of SNP effect
void gcta::blup_snp_geno()
{
    check_autosome();

    if(_mu.empty()) calcu_mu();

    int i=0, j=0, k=0, col_num=_varcmp_Py.cols();
    double x=0.0, fcount=0.0;

	// Calcuate A matrix
	cout<<"Calculating the BLUP solution to SNP effects ..."<<endl;
	vector<double> var_SNP(_include.size());
	eigenMatrix b_SNP=eigenMatrix::Zero(_include.size(), col_num); // variance of each SNP, 2pq
	for(j=0; j<_include.size(); j++){
        var_SNP[j]=_mu[_include[j]]*(1.0-0.5*_mu[_include[j]]);
        if(fabs(var_SNP[j])<1.0e-50) var_SNP[j]=0.0;
        else var_SNP[j]=1.0/var_SNP[j];
	}
	for(k=0; k<_include.size(); k++){
	    fcount=0.0;
        for(i=0; i<_keep.size(); i++){
            if(!_snp_1[_include[k]][i] || _snp_2[_include[k]][i]){
                if(_allele1[_include[k]]==_ref_A[_include[k]]) x=_snp_1[_include[k]][i]+_snp_2[_include[k]][i];
                else x=2.0-(_snp_1[_include[k]][i]+_snp_2[_include[k]][i]);
                x=(x-_mu[_include[k]]);
                for(j=0; j<col_num; j++) b_SNP(k,j)+=x*_varcmp_Py(i,j);
                fcount+=1.0;
            }
        }
        for(j=0; j<col_num; j++) b_SNP(k,j)=(b_SNP(k,j)*var_SNP[k]/fcount)*((double)_keep.size()/(double)_include.size());
	}
	output_blup_snp(b_SNP);
}

void gcta::blup_snp_dosage()
{
    check_autosome();

    if(_mu.empty()) calcu_mu();

    int i=0, j=0, k=0, col_num=_varcmp_Py.cols();

    // Substract each element by 2p
    for(i=0; i<_keep.size(); i++){
         for(j=0; j<_include.size(); j++) _geno_dose[i][_include[j]]-=_mu[_include[j]];
	}

	// Calcuate A matrix
	cout<<"Calculating the BLUP solution to SNP effects using imputed dosage scores ... "<<endl;
	vector<double> var_SNP(_include.size()); // variance of each SNP, 2pq
    eigenMatrix b_SNP=eigenMatrix::Zero(_include.size(), col_num); // variance of each SNP, 2pq
	for(j=0; j<_include.size(); j++){
        for(i=0; i<_keep.size(); i++) var_SNP[j]+=(double)_geno_dose[i][_include[j]]*_geno_dose[i][_include[j]];
        var_SNP[j]/=(double)(_keep.size()-1);
        if(fabs(var_SNP[j])<1.0e-50) var_SNP[j]=0.0;
        else var_SNP[j]=1.0/var_SNP[j];
	}
	for(k=0; k<_include.size(); k++){
	    for(i=0; i<_keep.size(); i++){
            for(j=0; j<col_num; j++) b_SNP(k,j)+=_geno_dose[i][_include[k]]*_varcmp_Py(i,j);
        }
        for(j=0; j<col_num; j++) b_SNP(k,j)=b_SNP(k,j)*var_SNP[k]/(double)_include.size();
	}
	output_blup_snp(b_SNP);
}

void gcta::output_blup_snp(eigenMatrix &b_SNP)
{
    string o_b_snp_file=_out+".snp.blp";
    ofstream o_b_snp(o_b_snp_file.c_str());
    if(!o_b_snp) throw("Error: can not open the file "+o_b_snp_file+" to write.");
    int i=0, j=0, col_num=b_SNP.cols();
	cout<<"Writing BLUP solutions of SNP effects for "<<_include.size()<<" SNPs to ["+o_b_snp_file+"]."<<endl;
    for(i=0; i<_include.size(); i++){
        o_b_snp<<_snp_name[_include[i]]<<"\t"<<_ref_A[_include[i]]<<"\t";
        for(j=0; j<col_num; j++) o_b_snp<<b_SNP(i,j)<<"\t";
        o_b_snp<<endl;
    }
    o_b_snp.close();
	cout<<"BLUP solutions of SNP effects for "<<_include.size()<<" SNPs have been saved in the file ["+o_b_snp_file+"]."<<endl;
}


