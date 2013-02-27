/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for estimating the genetic relationship matrix
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::check_autosome()
{
    for(int i=0; i<_include.size(); i++){
        if(_chr[_include[i]]>_autosome_num) throw("Error: this option is for the autosomal SNPs only.");
    }
}

void gcta::check_chrX()
{
   for(int i=0; i<_include.size(); i++){
        if(_chr[_include[i]]!=(_autosome_num+1)) throw("Error: this option is for SNPs on the X chromosome only.");
    }
}

void gcta::check_sex()
{
    for(int i=0; i<_keep.size(); i++){
        if(_sex[_keep[i]]!=1 && _sex[_keep[i]]!=2) throw("Error: Sex information of the individual \""+_fid[_keep[i]]+" "+_pid[_keep[i]]+"\" is missing.\nUse --update-sex option to update the sex information of the individuals.");
    }
}

// paa: proportion of ancestral alleles
void gcta::paa(string aa_file)
{
    check_autosome();
    if(_mu.empty()) calcu_mu();
	rm_high_ld();

    int i=0, j=0, k=0;

    // read ancestral alleles from a file
    ifstream i_aa(aa_file.c_str());
    if(!i_aa) throw("Error: can not open the file ["+aa_file+"] to read.");
    char cbuf='.';
    string str_buf;
	cout<<"Reading ancestral alleles of the SNPs from ["+aa_file+"]."<<endl;
    map<string, int>::iterator iter, End=_snp_name_map.end();
    vector<char> aa(_snp_num);
    for(i=0; i<_snp_num; i++) aa[i]='.';
    int icount=0;
	while(i_aa){
		i_aa>>str_buf;
		if(i_aa.eof()) break;
		iter=_snp_name_map.find(str_buf);
		i_aa>>cbuf;
		if(iter!=End && cbuf!='.'){
		    aa[iter->second]=cbuf;
		    icount++;
		}
	}
    i_aa.close();
	cout<<"Ancestral alleles for "<<icount<<" SNPs are included from ["+aa_file+"]."<<endl;

	cout<<"Calculating proportion of ancestral alleles ..."<<endl;
	double x=0.0;
	vector<double> hom_aa_rare(_keep.size()), hom_aa_comm(_keep.size()), hom_da_rare(_keep.size()), hom_da_comm(_keep.size()), het_aa_rare(_keep.size()), het_aa_comm(_keep.size()), nomiss(_keep.size());
	for(i=0; i<_keep.size(); i++){
 		for(k=0; k<_include.size(); k++){
 		    if(aa[_include[k]]=='.') continue;
            if(!_snp_1[_include[k]][i] || _snp_2[_include[k]][i]){
                x=_snp_1[_include[k]][i]+_snp_2[_include[k]][i];
                if(x<0.1){
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) hom_da_rare[i]+=1.0;
                        else hom_da_comm[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) hom_aa_rare[i]+=1.0;
                        else hom_aa_comm[i]+=1.0;
                    }
                }
                else if(x>1.9){
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) hom_aa_comm[i]+=1.0;
                        else hom_aa_rare[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) hom_da_comm[i]+=1.0;
                        else hom_da_rare[i]+=1.0;
                    }
                }
                else{
                    if(_ref_A[_include[k]]==aa[_include[k]]){
                        if(_mu[_include[k]]>1.0) het_aa_comm[i]+=1.0;
                        else het_aa_rare[i]+=1.0;
                    }
                    else{
                        if(_mu[_include[k]]>1.0) het_aa_rare[i]+=1.0;
                        else het_aa_comm[i]+=1.0;
                    }
                }
                nomiss[i]+=1.0;
            }
        }
        hom_aa_rare[i]/=nomiss[i];
        hom_aa_comm[i]/=nomiss[i];
        hom_da_rare[i]/=nomiss[i];
        hom_da_comm[i]/=nomiss[i];
        het_aa_rare[i]/=nomiss[i];
        het_aa_comm[i]/=nomiss[i];
        cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}

    // Save matrix A in binary file
	string paa_file=_out+".paa";
	ofstream o_paa(paa_file.c_str());
	if(!o_paa) throw("Error: can not open the file ["+paa_file+"] to write.");
	o_paa<<"FID\tIID\tNOMISS\tHOM_AA_RARE\tHOM_AA_COMM\tHOM_DA_RARE\tHOM_DA_COMM\tHET_AA_RARE\tHET_AA_COMM"<<endl;
	for(i=0; i<_keep.size(); i++) o_paa<<_fid[i]<<"\t"<<_pid[i]<<"\t"<<nomiss[i]<<"\t"<<hom_aa_rare[i]<<"\t"<<hom_aa_comm[i]<<"\t"<<hom_da_rare[i]<<"\t"<<hom_da_comm[i]<<"\t"<<het_aa_rare[i]<<"\t"<<het_aa_comm[i]<<endl;
	o_paa.close();
	cout<<"Proportion of ancestral alleles has been saved in file ["+paa_file+"]."<<endl;
}

// inbreeding coefficient
void gcta::ibc(bool ibc_all)
{
    check_autosome();
    if(_mu.empty()) calcu_mu();
	rm_high_ld();

    int i=0, j=0, k=0;

	// Calcuate A matrix
	cout<<"Calculating the inbreeding coefficients ..."<<endl;
	double x=0.0, sum_w=0.0, sum_h=0.0, Fhat_buf=0.0;
	vector<double> h(_include.size()), h_i(_include.size()), w(_include.size()), p(_include.size()), p_q(_include.size()), q_p(_include.size()); // variance of each SNP, 2pq
	for(j=0; j<_include.size(); j++){
	    p[j]=0.5*_mu[_include[j]];
        h[j]=2.0*p[j]*(1.0-p[j]);
        p_q[j]=p[j]/(1.0-p[j]);
        if(fabs(p_q[j])<1.0e-50) q_p[j]=0.0;
        else q_p[j]=1.0/p_q[j];
        w[j]=h[j]/(1.0-h[j]);
        if(fabs(h[j])<1.0e-50) h_i[j]=0.0;
        else h_i[j]=1.0/h[j];
	}
	vector<double> Fhat1(_keep.size()), Fhat1_w(_keep.size()), Fhat2(_keep.size()), Fhat2_w(_keep.size()), Fhat3(_keep.size()), Fhat4(_keep.size()), Fhat5(_keep.size()), Fhat6(_keep.size()), Fhat7(_keep.size());
	vector<double> rare_hom(_keep.size()), comm_hom(_keep.size()), nomiss(_keep.size());
	for(i=0; i<_keep.size(); i++){
        sum_w=0.0;
        sum_h=0.0;
		for(k=0; k<_include.size(); k++){
            if(!_snp_1[_include[k]][i] || _snp_2[_include[k]][i]){
                x=_snp_1[_include[k]][i]+_snp_2[_include[k]][i];
                if(_allele2[_include[k]]==_ref_A[_include[k]]) x=2.0-x;
                Fhat_buf=(x-_mu[_include[k]])*(x-_mu[_include[k]]);
                Fhat4[i]+=Fhat_buf;
                Fhat_buf*=h_i[k];
                Fhat1[i]+=Fhat_buf;
                Fhat1_w[i]+=Fhat_buf*w[k];
                Fhat_buf=h[k]-x*(2.0-x);
                Fhat6[i]+=Fhat_buf;
                Fhat_buf*=h_i[k];
                Fhat2[i]+=Fhat_buf;
                Fhat2_w[i]+=Fhat_buf*w[k];
                Fhat_buf=(x*(x-1.0-_mu[_include[k]])+_mu[_include[k]]*p[k]);
                Fhat3[i]+=Fhat_buf*h_i[k];
                Fhat5[i]+=Fhat_buf;
                sum_w+=w[k];
                sum_h+=h[k];
                if(x<0.1) Fhat7[i]+=p_q[k];
                else if(x>1.9) Fhat7[i]+=q_p[k];
                // Count the number of rare and common homozygotes
                if(x<0.1){
                    if(p[k]>0.5) rare_hom[i]+=1.0;
                    else comm_hom[i]+=1.0;
                }
                else if(x>1.9){
                    if(p[k]<0.5) rare_hom[i]+=1.0;
                    else comm_hom[i]+=1.0;
                }
                nomiss[i]+=1.0;
            }
        }
        Fhat1[i]/=nomiss[i];
        Fhat2[i]/=nomiss[i];
        Fhat3[i]/=nomiss[i];
        Fhat1_w[i]/=sum_w;
        Fhat2_w[i]/=sum_w;
        Fhat4[i]/=sum_h;
        Fhat5[i]/=sum_h;
        Fhat6[i]/=sum_h;
        rare_hom[i]/=nomiss[i];
        comm_hom[i]/=nomiss[i];
        Fhat7[i]/=nomiss[i];
        cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}

    // Save matrix A in binary file
	string ibc_file=_out+".ibc";
	ofstream o_ibc(ibc_file.c_str());
	if(!o_ibc) throw("Error: can not open the file ["+ibc_file+"] to write.");
	if(ibc_all){
        o_ibc<<"FID\tIID\tNOMISS\tP_RARE_HOM\tP_COMM_HOM\tFhat1\tFhat1_w\tFhat2\tFhat2_w\tFhat3\tFhat4\tFhat5\tFhat6\tFhat7"<<endl;
        for(i=0; i<_keep.size(); i++) o_ibc<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t"<<nomiss[i]<<"\t"<<rare_hom[i]<<"\t"<<comm_hom[i]<<"\t"<<Fhat1[i]-1.0<<"\t"<<Fhat1_w[i]-1.0<<"\t"<<Fhat2[i]<<"\t"<<Fhat2_w[i]<<"\t"<<Fhat3[i]<<"\t"<<Fhat4[i]-1.0<<"\t"<<Fhat5[i]<<"\t"<<Fhat6[i]<<"\t"<<Fhat7[i]<<endl;
	}
	else{
        o_ibc<<"FID\tIID\tNOMISS\tFhat1\tFhat2\tFhat3"<<endl;
        for(i=0; i<_keep.size(); i++) o_ibc<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t"<<nomiss[i]<<"\t"<<Fhat1[i]-1.0<<"\t"<<Fhat2[i]<<"\t"<<Fhat3[i]<<endl;
	}
	o_ibc.close();
	cout<<"Inbreeding coefficients have been saved in the file ["+ibc_file+"]."<<endl;
}

void gcta::make_grm(bool grm_xchr_flag, bool output_bin, int grm_mtd)
{
    if(grm_xchr_flag) check_chrX();
    else check_autosome();
	rm_high_ld();

	int i=0, j=0, k=0;
    vector< vector<float> > X;
    make_XMat(X, false);
    vector<double> sd_SNP_i, d(_include.size());
    std_XMat(X, sd_SNP_i, grm_xchr_flag, !grm_mtd);
    if(grm_mtd==1){
        for(j=0; j<_include.size(); j++){
            for(i=0, k=0; i<_keep.size(); i++){
                if(X[i][j]<1e5){ d[j]+=X[i][j]*X[i][j]; k++; }
            }
            if(k>0) d[j]/=(double)k;
        }
    }

	// Calcuate A matrix
	if(_dosage_flag){
        if(grm_xchr_flag) cout<<"\nCalculating the genetic relationship matrix for the X chromosome using imputed dosage data ... "<<endl;
	    else cout<<"\nCalculating the genetic relationship matrix using imputed dosage data ... "<<endl;
	}
	else{
        if(grm_xchr_flag) cout<<"\nCalculating the genetic relationship matrix for the X chromosome ... ((NOTE: default speed-optimized mode, may use huge RAM))"<<endl;
        else cout<<"\nCalculating the genetic relationship matrix ... (NOTE: default speed-optimized mode, may use huge RAM)"<<endl;
	}
	double denom=0.0;
	vector< vector<double> > A(_keep.size());
	vector< vector<int> > A_N(_keep.size());
	for(i=0; i<_keep.size(); i++){
	    A[i].resize(_keep.size());
	    A_N[i].resize(_keep.size());
	}
	for(i=0; i<_keep.size(); i++){
	    if(grm_mtd==1){
            for(k=0, denom=0.0; k<_include.size(); k++){
                if(X[i][k]<1e5){
                    A[i][i]+=X[i][k]*X[i][k];
                    denom+=d[k];
                    A_N[i][i]++;
                }
            }
            if(A_N[i][i]>0) A[i][i]/=denom;
	    }
	    else{
            for(k=0; k<_include.size(); k++){
                if(X[i][k]<1e5){
                    A[i][i]+=X[i][k]*X[i][k];
                    A_N[i][i]++;
                }
            }
            if(A_N[i][i]>0) A[i][i]/=(double)(A_N[i][i]);
	    }
	    if(A_N[i][i]==0) A[i][i]=1.0;
		for(j=0; j<i; j++){
		    if(grm_mtd==1){
                for(k=0, denom=0.0; k<_include.size(); k++){
                    if(X[i][k]<1e5 && X[j][k]<1e5){
                        A[i][j]+=X[i][k]*X[j][k];
                        denom+=d[k];
                        A_N[i][j]++;
                    }
                }
                if(A_N[i][j]>0.0) A[i][j]/=denom;
		    }
		    else{
                for(k=0; k<_include.size(); k++){
                    if(X[i][k]<1e5 && X[j][k]<1e5){
                        A[i][j]+=X[i][k]*X[j][k];
                        A_N[i][j]++;
                    }
                }
                if(A_N[i][j]>0.0) A[i][j]/=(double)A_N[i][j];
		    }
		}
		cout<<i+1<<" of "<<_keep.size()<<" individuals.\r";
	}
	cout<<endl<<endl;
	output_grm_vec(A, A_N, output_bin);
}

void gcta::output_grm_vec(vector< vector<double> > &A, vector< vector<int> > &A_N, bool output_grm_bin)
{
    int i=0, j=0;
    string grm_file;
    if(output_grm_bin){
        // Save matrix A in binary file
        grm_file=_out+".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out|ios::binary);
        if(!A_Bin) throw("Error: can not open the file ["+grm_file+"] to write.");
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<i; j++) A_Bin.write((char *) &(A[i][j]), sizeof(double));
            for(j=i; j<_keep.size(); j++) A_Bin.write((char *) &(A[j][i]), sizeof(double));
        }
        A_Bin.close();
        cout<<"GRM of "<<_keep.size()<<" individuals has been saved in the file ["+grm_file+"] (in binary format)."<<endl;
    }
    else{
        // Save A matrix in txt format
        grm_file=_out+".grm.gz";
        gzofstream zoutf;
        zoutf.open( grm_file.c_str() );
        if(!zoutf.is_open()) throw("Error: can not open the file ["+grm_file+"] to write.");
        cout<<"Saving the genetic relationship matrix to the file ["+grm_file+"] (in compressed text format)."<<endl;
        zoutf.setf(ios::scientific);
        zoutf.precision(6);
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) zoutf<<i+1<<'\t'<<j+1<<'\t'<<A_N[i][j]<<'\t'<<A[i][j]<<endl;
        }
        zoutf.close();
        cout<<"The genetic relationship matrix has been saved in the file ["+grm_file+"] (in compressed text format)."<<endl;
    }

    string famfile=_out+".grm.id";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the file ["+famfile+"] to write.");
	for(i=0; i<_keep.size(); i++) Fam<<_fid[_keep[i]]+"\t"+_pid[_keep[i]]<<endl;
	Fam.close();
	cout<<"IDs for the GRM file ["+grm_file+"] have been saved in the file ["+famfile+"]."<<endl;
}

void gcta::output_grm_MatrixXf(bool output_grm_bin)
{
    int i=0, j=0;
    string grm_file;
    if(output_grm_bin){
        // Save matrix A in binary file
        grm_file=_out+".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out|ios::binary);
        if(!A_Bin) throw("Error: can not open the file ["+grm_file+"] to write.");
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<i; j++) A_Bin.write((char *) &(_grm(i,j)), sizeof(double));
            for(j=i; j<_keep.size(); j++) A_Bin.write((char *) &(_grm(j,i)), sizeof(double));
        }
        A_Bin.close();
        cout<<"GRM of "<<_keep.size()<<" individuals has been saved in the file ["+grm_file+"] (in binary format)."<<endl;
    }
    else{
        // Save A matrix in txt format
        grm_file=_out+".grm.gz";
        gzofstream zoutf;
        zoutf.open( grm_file.c_str() );
        if(!zoutf.is_open()) throw("Error: can not open the file ["+grm_file+"] to write.");
        cout<<"Saving the genetic relationship matrix to the file ["+grm_file+"] (in compressed text format)."<<endl;
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) zoutf<<setprecision(0)<<i+1<<'\t'<<j+1<<'\t'<<(int)_grm_N(i,j)<<'\t'<<setprecision(15)<<_grm(i,j)<<endl;
        }
        zoutf.close();
        cout<<"The genetic relationship matrix has been saved in the file ["+grm_file+"] (in compressed text format)."<<endl;
    }

    string famfile=_out+".grm.id";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the file ["+famfile+"] to write.");
	for(i=0; i<_keep.size(); i++) Fam<<_fid[_keep[i]]+"\t"+_pid[_keep[i]]<<endl;
	Fam.close();
	cout<<"IDs for the GRM file ["+grm_file+"] have been saved in the file ["+famfile+"]."<<endl;
}

void gcta::read_grm_gz(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
    // read GRM IDs
    string grm_id_file=grm_file+".grm.id";
    if(out_id_log) cout<<"Reading IDs of the genetic relationship matrix (GRM) from ["+grm_id_file+"]."<<endl;
    ifstream i_grm_id(grm_id_file.c_str());
	if(!i_grm_id) throw("Error: can not open the file ["+grm_id_file+"] to read.");
	string str_buf, id_buf;
	vector<string> fid, pid;
	grm_id.clear();
	while(i_grm_id){
	    i_grm_id>>str_buf;
	    if(i_grm_id.eof()) break;
	    fid.push_back(str_buf);
		id_buf=str_buf+":";
	    i_grm_id>>str_buf;
	    pid.push_back(str_buf);
		id_buf+=str_buf;
		grm_id.push_back(id_buf);
		getline(i_grm_id, str_buf);
	}
	i_grm_id.close();
	int n=grm_id.size();
    if(out_id_log) cout<<n<<" IDs read from ["+grm_id_file+"]."<<endl;

    if(_id_map.empty()){
        _fid=fid;
        _pid=pid;
        _indi_num=_fid.size();
        _sex.resize(_fid.size());
        init_keep();
    }

    if(read_id_only) return;

    string grm_gzfile=grm_file+".grm.gz";
    const int MAX_LINE_LENGTH = 1000;
    char buf[MAX_LINE_LENGTH];
    gzifstream zinf;
    zinf.open(grm_gzfile.c_str());
    if(!zinf.is_open()) throw("Error: can not open the file ["+grm_gzfile+"] to read.");

    int indx1=0, indx2=0, nline=0;
    double grm_buf=0.0, grm_N_buf;
    string errmsg="Error: failed to read ["+grm_gzfile+"]. The format of the GRM file has been changed?\nError occurs in line:\n";
    cout<<"Reading the GRM from ["+grm_gzfile+"]."<<endl;
    //zinf.getline(buf, MAX_LINE_LENGTH, '\n'); // skip the header
    _grm.resize(n,n);
    _grm_N.resize(n,n);
    while(1){
        zinf.getline(buf, MAX_LINE_LENGTH, '\n');
        if(zinf.fail() || !zinf.good()) break;
        stringstream ss(buf);
        if(!(ss>>indx1)) throw(errmsg+buf);
        if(!(ss>>indx2)) throw(errmsg+buf);
        if(!(ss>>grm_N_buf)) throw(errmsg+buf);
        if(!(ss>>grm_buf)) throw(errmsg+buf);
		if(indx1 < indx2 || indx1>n || indx2>n) throw(errmsg+buf);
		if(grm_N_buf==0) cout<<"Warning: "<<buf<<endl;
		_grm_N(indx1-1,indx2-1)=grm_N_buf;
		_grm(indx1-1,indx2-1)=grm_buf;
		nline++;
        if(ss>>str_buf) throw(errmsg+buf);
    }
    zinf.close();
    if(nline!=(int)(n*(n+1)*0.5)) throw("Error: incorrect number of lines in the grm file. *.grm.gz file and *.grm.id file are mismatched?");
    cout<<"Pairwise genetic relationships between "<<n<<" individuals are included from ["+grm_gzfile+"]."<<endl;
}

void gcta::rm_cor_indi(double grm_cutoff)
{
    cout<<"Pruning the GRM with a cutoff of "<<grm_cutoff<<" ..."<<endl;

    int i=0, j=0, i_buf=0;

    vector<int> rm_grm_ID1, rm_grm_ID2;
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<i; j++){
            if(_grm(_keep[i],_keep[j])>grm_cutoff){
                rm_grm_ID1.push_back(_keep[i]);
                rm_grm_ID2.push_back(_keep[j]);
            }
        }
    }
    vector<int> rm_uni_ID(rm_grm_ID1);
	rm_uni_ID.insert(rm_uni_ID.end(), rm_grm_ID2.begin(), rm_grm_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    map<int,int> rm_uni_ID_count;
    for(i=0; i<rm_uni_ID.size(); i++){
        i_buf=count(rm_grm_ID1.begin(), rm_grm_ID1.end(), rm_uni_ID[i])+count(rm_grm_ID2.begin(), rm_grm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(pair<int,int>(rm_uni_ID[i], i_buf));
    }
    map<int,int>::iterator iter1, iter2;
    for(i=0; i<rm_grm_ID1.size(); i++){
        iter1=rm_uni_ID_count.find(rm_grm_ID1[i]);
        iter2=rm_uni_ID_count.find(rm_grm_ID2[i]);
        if(iter1->second < iter2->second){
            i_buf=rm_grm_ID1[i];
            rm_grm_ID1[i]=rm_grm_ID2[i];
            rm_grm_ID2[i]=i_buf;
        }
    }
    stable_sort(rm_grm_ID1.begin(), rm_grm_ID1.end());
    rm_grm_ID1.erase(unique(rm_grm_ID1.begin(), rm_grm_ID1.end()), rm_grm_ID1.end());
    vector<string> removed_ID;
    for(i=0; i<rm_grm_ID1.size(); i++) removed_ID.push_back(_fid[rm_grm_ID1[i]]+":"+_pid[rm_grm_ID1[i]]);

    // update _keep and _id_map
    update_id_map_rm(removed_ID, _id_map, _keep);

    cout<<"After pruning the GRM, there are "<<_keep.size()<<" individuals ("<<removed_ID.size()<<" individuals removed)."<<endl;
}

void gcta::adj_grm(double adj_grm_fac)
{
    cout<<"Adjusting the GRM for sampling errors ..."<<endl;
    int i=0, j=0, n=_keep.size();
    double off_mean=0.0, diag_mean=0.0, off_var=0.0, diag_var=0.0, d_buf=0.0;
    for(i=0; i<n; i++){
        diag_mean+=_grm(_keep[i],_keep[i]);
        for(j=0; j<i; j++) off_mean+=_grm(_keep[i],_keep[j]);
    }
    diag_mean/=n;
    off_mean/=0.5*n*(n-1.0);
    for(i=0; i<n; i++){
        d_buf=_grm(_keep[i],_keep[i])-diag_mean;
        diag_var+=d_buf*d_buf;
        for(j=0; j<i; j++){
            d_buf=_grm(_keep[i],_keep[j])-off_mean;
            off_var+=d_buf*d_buf;
        }
    }
    diag_var/=n-1.0;
    off_var/=0.5*n*(n-1.0)-1.0;
    for(i=0; i<_keep.size(); i++){
        d_buf=1.0-(adj_grm_fac+1.0/_grm_N(_keep[i],_keep[i]))/diag_var;
        if(_grm(_keep[i],_keep[i])>0) _grm(_keep[i],_keep[i])=1.0+d_buf*(_grm(_keep[i],_keep[i])-1.0);
        for(j=0; j<i; j++){
            if(_grm_N(_keep[i],_keep[j])>0) _grm(_keep[i],_keep[j])*=1.0-(adj_grm_fac+1.0/_grm_N(_keep[i],_keep[j]))/off_var;
        }
    }
}

void gcta::dc(int dosage_compen)
{
    cout<<"Parameterizing the GRM under the assumption of ";
    if(dosage_compen==1) cout<<"full dosage compensation ..."<<endl;
    else if(dosage_compen==0) cout<<"no dosage compensation ..."<<endl;

    int i=0, j=0, i_buf=0;
    double c1=1.0, c2=1.0;
    if(dosage_compen==1){ c1=2.0; c2=sqrt(2.0); } // full dosage compensation
    else if(dosage_compen==0){ c1=0.5; c2=sqrt(0.5); } // on dosage compensation
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<=i; j++){
            i_buf=_sex[_keep[i]]*_sex[_keep[j]];
            if(i_buf==1) _grm(i,j)*=c1;
            else if(i_buf==2) _grm(i,j)*=c2;
        }
    }
}

void gcta::manipulate_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag)
{
    int i=0, j=0;

    vector<string> grm_id;
    if(merge_grm_flag) merge_grm(grm_file);
    else read_grm_gz(grm_file, grm_id);

    if(!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if(!remove_indi_file.empty()) remove_indi(remove_indi_file);
    if(grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
    if(!sex_file.empty()) update_sex(sex_file);
    if(adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
    if(dosage_compen>-1) dc(dosage_compen);
    if(grm_cutoff>-1.0 || !keep_indi_file.empty() || !remove_indi_file.empty()){
        eigenMatrix tmp(_grm);
        _grm.resize(_keep.size(),_keep.size());
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) _grm(i,j)=tmp(_keep[i],_keep[j]);
        }
        tmp=_grm_N;
        _grm_N.resize(_keep.size(),_keep.size());
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) _grm_N(i,j)=tmp(_keep[i],_keep[j]);
        }
    }
}

void gcta::save_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool output_grm_bin)
{
    if(dosage_compen>-1) check_sex();
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, sex_file, grm_cutoff, adj_grm_fac, dosage_compen, merge_grm_flag);
    output_grm_MatrixXf(output_grm_bin);
}

void gcta::pca(string grm_file, string keep_indi_file, string remove_indi_file, double grm_cutoff, bool merge_grm_flag, int out_pc_num)
{
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, "", grm_cutoff, -2.0, -2, merge_grm_flag);
    _grm_N.resize(1,1);
    int i=0, j=0;
    cout<<"\nPerforming principal component analysis ..."<<endl;
    SelfAdjointEigenSolver<eigenMatrix> eigensolver(_grm);
    eigenMatrix evec=(eigensolver.eigenvectors());
    eigenVector eval=eigensolver.eigenvalues();

    string eval_file=_out+".eigenval";
    ofstream o_eval(eval_file.c_str());
    if(!o_eval) throw("Error: can not open the file ["+eval_file+"] to read.");
    for(i=_keep.size()-1; i>=0; i--) o_eval<<eval(i)<<endl;
    o_eval.close();
    cout<<"Eigenvalues of "<<_keep.size()<<" individuals have been saved in ["+eval_file+"]."<<endl;
    string evec_file=_out+".eigenvec";
    ofstream o_evec(evec_file.c_str());
    if(!o_evec) throw("Error: can not open the file ["+evec_file+"] to read.");
    if(out_pc_num>_keep.size()) out_pc_num=_keep.size();
    for(i=0; i<_keep.size(); i++){
        o_evec<<_fid[_keep[i]]<<" "<<_pid[_keep[i]];
        for(j=_keep.size()-1; j>=(_keep.size()-out_pc_num); j--) o_evec<<" "<<evec(i,j);
        o_evec<<endl;
    }
    o_evec.close();
    cout<<"The first "<<out_pc_num<<" eigenvectors of "<<_keep.size()<<" individuals have been saved in ["+evec_file+"]."<<endl;
}

void gcta::merge_grm(string merge_grm_file)
{
    vector<string> grm_files, grm_id;
    read_grm_filenames(merge_grm_file, grm_files);

    int f=0, i=0, j=0;
    for(f=0; f<grm_files.size(); f++){
        read_grm_gz(grm_files[f], grm_id, false, true);
        update_id_map_kp(grm_id, _id_map, _keep);
    }
    vector<string> uni_id;
	for(i=0; i<_keep.size(); i++) uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
	_n=uni_id.size();
	if(_n==0) throw("Error: no individual is in common in the GRM files.");
	else cout<<_n<<" individuals in common in the GRM files."<<endl;

    vector<int> kp;
    eigenMatrix grm=eigenMatrix::Zero(_n, _n);
    eigenMatrix grm_N=eigenMatrix::Zero(_n, _n);
    for(f=0; f<grm_files.size(); f++){
        cout<<"Reading the GRM from the "<<f+1<<"th file ..."<<endl;
        read_grm_gz(grm_files[f], grm_id);
        StrFunc::match(uni_id, grm_id, kp);
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++){
                if(kp[i]>=kp[j]){
                    grm(i,j)+=_grm(kp[i],kp[j])*_grm_N(kp[i],kp[j]);
                    grm_N(i,j)+=_grm_N(kp[i],kp[j]);
                }
                else{
                    grm(i,j)+=_grm(kp[j],kp[i])*_grm_N(kp[j],kp[i]);
                    grm_N(i,j)+=_grm_N(kp[j],kp[i]);
                }
            }
        }
    }
    for(i=0; i<_n; i++){
        for(j=0; j<=i; j++){
            if(grm_N(i,j)==0) _grm(i,j)=0;
            else _grm(i,j)=grm(i,j)/grm_N(i,j);
            _grm_N(i,j)=grm_N(i,j);
        }
    }
    grm.resize(1,1);
    grm_N.resize(1,1);
    cout<<"\n"<<grm_files.size()<<" GRMs have been merged together."<<endl;
}

void gcta::read_grm_filenames(string merge_grm_file, vector<string> &grm_files, bool out_log)
{
    ifstream merge_grm(merge_grm_file.c_str());
    if(!merge_grm) throw("Error: can not open the file ["+merge_grm_file+"] to read.");
    string str_buf;
    grm_files.clear();
    vector<string> vs_buf;
    while(getline(merge_grm, str_buf)){
        if(!str_buf.empty()){
            if(StrFunc::split_string(str_buf, vs_buf)==1) grm_files.push_back(vs_buf[0]);
        }
    }
    if(out_log) cout<<"There are "<<grm_files.size()<<" GRM file names specified in ["+merge_grm_file+"]."<<endl;
    if(grm_files.size()>1000) throw("Error: too many GRM file names specified in ["+merge_grm_file+"]. Maximum is 1000.");
    if(grm_files.size()<1) throw("Error: no GRM file name is found in ["+merge_grm_file+"].");
}
