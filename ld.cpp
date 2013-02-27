/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for estimating the LD structure
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_LD_target_SNPs(string snplistfile)
{
    // Read snplist file
    _ld_target_snp.clear();
	read_snplist(snplistfile, _ld_target_snp);

	int i=0, prev_target_snp_num=_ld_target_snp.size(), prev_include_size=_include.size();
	map<string, int> snp_map_buf;
	map<string, int>::iterator iter;
	for(i=0; i<_snp_num; i++) snp_map_buf.insert(pair<string, int>(_snp_name[i], i));
	for(i=_ld_target_snp.size()-1; i>=0; i--){
	    iter=snp_map_buf.find(_ld_target_snp[i]);
	    if(iter!=snp_map_buf.end()) _snp_name_map.insert(*iter);
	    else _ld_target_snp.erase(_ld_target_snp.begin()+i);
	}
    _include.clear();
    for(iter=_snp_name_map.begin(); iter!=_snp_name_map.end(); iter++) _include.push_back(iter->second);
    stable_sort(_include.begin(), _include.end());
	if(_ld_target_snp.size()==0) throw("Error: no target SNPs are retained to estimate the LD structure.");
	else cout<<prev_target_snp_num<<" target SNPs read from ["+snplistfile+"], "<<_ld_target_snp.size()<<" of which exist in the data."<<endl;
}

void gcta::LD_Blocks(int stp, double wind_size, double alpha, bool IncldQ, bool save_ram)
{
	int i=0, j=0;
	if(_mu.empty()) calcu_mu();

    // Read snplist file
 	vector<int> smpl_buf, smpl;
	vector<string> uni_snp;
	for(i=0; i<_include.size(); i++) uni_snp.push_back(_snp_name[_include[i]]);
	StrFunc::match(_ld_target_snp, uni_snp, smpl_buf);
	for(i=0; i<smpl_buf.size(); i++){
	    if(smpl_buf[i]>-1) smpl.push_back(smpl_buf[i]);
	}
    int SNP_SmplNum=smpl.size(); // smpl is the position of _include

    // Calculate LD structure
    cout<<"Estimating LD structure..."<<endl;
    vector<int> K(SNP_SmplNum);
	vector<double> r2(SNP_SmplNum), md_r2(SNP_SmplNum), max_r2(SNP_SmplNum), dL(SNP_SmplNum), dR(SNP_SmplNum);
	vector<string> max_r2_snp(SNP_SmplNum);
	vector< vector<double> > r(SNP_SmplNum);
	vector<string> L_SNP(SNP_SmplNum), R_SNP(SNP_SmplNum);
	vector< vector<string> > snp_ls(SNP_SmplNum);
	EstLD(smpl, wind_size, snp_ls, r, r2, md_r2, max_r2, max_r2_snp, dL, dR, K, L_SNP, R_SNP, alpha, IncldQ);

	// Save result
	string SavFileName=_out+".rsq.ld";
	ofstream SavFile(SavFileName.c_str());
	if(!SavFile) throw("Error: can not open the file ["+SavFileName+"] to save result!");
	SavFile<<"target_SNP\tfreq\tL_region\tR_region\tL_snp\tR_snp\tnSNPs\tmean_rsq\tmedian_rsq\tmax_rsq\tmax_rsq_snp"<<endl;
	for(i=0; i<SNP_SmplNum; i++) SavFile<<_snp_name[_include[smpl[i]]]<<"\t"<<0.5*_mu[_include[smpl[i]]]<<"\t"<<dL[i]<<"\t"<<dR[i]<<"\t"<<L_SNP[i]<<"\t"<<R_SNP[i]<<"\t"<<K[i]<<"\t"<<r2[i]<<"\t"<<md_r2[i]<<"\t"<<max_r2[i]<<"\t"<<max_r2_snp[i]<<endl;
	SavFile.close();
	SavFileName=_out+".r.ld";
	SavFile.open(SavFileName.c_str());
	if(!SavFile) throw("Error: can not open the file ["+SavFileName+"] to save result.");
	for(i=0; i<SNP_SmplNum; i++){
        for(j=0; j<r[i].size(); j++) SavFile<<r[i][j]<<" ";
        SavFile<<endl;
    }
	SavFile.close();
	SavFileName=_out+".snp.ld";
	SavFile.open(SavFileName.c_str());
	if(!SavFile) throw("Can not open the file ["+SavFileName+"] to save result.");
	for(i=0; i<SNP_SmplNum; i++){
        for(j=0; j<snp_ls[i].size(); j++) SavFile<<snp_ls[i][j]<<" ";
        SavFile<<endl;
    }
	SavFile.close();
	cout<<"Results have been saved in ["+_out+".rsq.ld]"+", ["+_out+".r.ld]"+" and ["+_out+".snp.ld].\n"<<endl;
}

void gcta::EstLD(vector<int> &smpl, double wind_size, vector< vector<string> > &snp, vector< vector<double> > &r, vector<double> &r2, vector<double> &md_r2, vector<double> &max_r2, vector<string> &max_r2_snp, vector<double> &dL, vector<double> &dR, vector<int> &K, vector<string> &L_SNP, vector<string> &R_SNP, double alpha, bool IncldQ)
{
	int i=0, j=0, L=0, R=0, maxL=0, maxR=0, i_buf=0;
    map<int, int> smpl_snp_map;
    for(i=0; i<smpl.size(); i++) smpl_snp_map.insert(pair<int,int>(_include[smpl[i]], i));

    cout<<"Parameters used to search SNPs in LD with the given SNPs: window size="<<(int)(wind_size*0.001)<<"Kb, significant level="<<alpha<<endl;
	vector<double> rst, y, x;
    for(i=0; i<smpl.size(); i++){
        vector<int> buf;
        vector<double> r_buf, rsq_buf;
        maxL=maxR=L=R=smpl[i];
		makex(L, y);
		if(IncldQ){ buf.push_back(L); rsq_buf.push_back(1.0); r_buf.push_back(1.0); }
        while(1){
            if(R==_include.size()-1) break;
			if(_chr[_include[R]]!=_chr[_include[R+1]]) break;
			if(_bp[_include[R+1]]-_bp[_include[smpl[i]]]>wind_size) break;
			R++;
			if(smpl_snp_map.find(_include[R])!=smpl_snp_map.end()) continue;
			makex(R, x);
			reg(y, x, rst);
			if(rst[2]<alpha){ maxR=R; buf.push_back(R); rsq_buf.push_back(rst[3]); r_buf.push_back(rst[4]); }
        }
        while(1){
            if(L==0) break;
			if(_chr[_include[L]]!=_chr[_include[L-1]]) break;
			if(_bp[_include[smpl[i]]]-_bp[_include[L-1]]>wind_size) break;
            L--;
            if(smpl_snp_map.find(_include[L])!=smpl_snp_map.end()) continue;
			makex(L, x);
			reg(y, x, rst);
			if(rst[2]<alpha){ maxL=L; buf.insert(buf.begin(), L); rsq_buf.insert(rsq_buf.begin(), rst[3]); r_buf.insert(r_buf.begin(), rst[4]); }
        }
        if(buf.size()==0){
            K[i]=0; dL[i]=0; dR[i]=0; L_SNP[i]="NA"; R_SNP[i]="NA";
            r[i].push_back(0.0); snp[i].push_back("NA");
            r2[i]=0.0; md_r2[i]=0.0; max_r2[i]=0.0; max_r2_snp[i]="NA";
        }
        else{
            K[i]=buf.size();
            dL[i]=_bp[_include[smpl[i]]]-_bp[_include[maxL]];
            dR[i]=_bp[_include[maxR]]-_bp[_include[smpl[i]]];
            L_SNP[i]=_snp_name[_include[maxL]];
            R_SNP[i]=_snp_name[_include[maxR]];
            for(j=0; j<K[i]; j++){
                r[i].push_back(r_buf[j]);
                snp[i].push_back(_snp_name[_include[buf[j]]]);
            }
            r2[i]=CommFunc::mean(rsq_buf);
            md_r2[i]=CommFunc::median(rsq_buf);
            i_buf=max_element(rsq_buf.begin(), rsq_buf.end())-rsq_buf.begin();
            max_r2[i]=rsq_buf[i_buf];
            max_r2_snp[i]=snp[i][i_buf];
        }
        cout<<i+1<<" of "<<smpl.size()<<" target SNPs.\r";
    }
}

void gcta::reg(vector<double> &y, vector<double> &x, vector<double> &rst)
{
    int N=x.size();
    if(N!=y.size() || N<1) throw("Error: The lengths of x and y do not match.");

    int i=0;
    double d_buf=0.0, y_mu=0.0, x_mu=0.0, x_var=0.0, y_var=0.0, cov=0.0;
    for(i=0; i<N; i++){ x_mu+=x[i]; y_mu+=y[i]; }
    x_mu/=(double)N;
    y_mu/=(double)N;
    for(i=0; i<N; i++){
        d_buf=(x[i]-x_mu);
        x_var+=d_buf*d_buf;
        d_buf=(y[i]-y_mu);
        y_var+=d_buf*d_buf;
    }
    x_var/=(double)(N-1.0);
    y_var/=(double)(N-1.0);
    for(i=0; i<N; i++) cov+=(x[i]-x_mu)*(y[i]-y_mu);
    cov/=(double)(N-1);

    double a=0.0, b=0.0, sse=0.0, b_se=0.0, p=0.0, rsq=0.0, r=0.0;
    if(x_var>0.0) b=cov/x_var;
    a=y_mu-b*x_mu;
    for(i=0; i<N; i++){
        d_buf=y[i]-a-b*x[i];
        sse+=d_buf*d_buf;
    }
    if(x_var>0.0) b_se=sqrt(sse/x_var/(N-1.0)/(N-2.0));
    if(x_var>0.0 && y_var>0.0){
        r=cov/sqrt(y_var*x_var);
        rsq=r*r;
    }
    double t=0.0;
    if(b_se>0.0) t=fabs(b/b_se);
    p=StatFunc::t_prob(N-2.0, t, true);
    rst.clear();
    rst.push_back(b); rst.push_back(b_se); rst.push_back(p); rst.push_back(rsq); rst.push_back(r);
}
