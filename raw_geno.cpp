/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for reading the raw genotype data
 *
 * 2010 by Jian Yang <jian.yang@qimr.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::read_IRG_fnames(string snp_info_file, string fname_file, double GC_cutoff)
{
    // read SNP summary information
    cout<<"Reading summary information of the SNPs ..."<<endl;
    ifstream i_snp_info(snp_info_file.c_str());
    if(!i_snp_info) throw("Error: can not open the file ["+snp_info_file+"] to read.");
    int i=0,j=0;
    string str_buf;
    vector<string> vs_buf;
    getline(i_snp_info, str_buf);
    while(getline(i_snp_info, str_buf)){
        StrFunc::split_string(str_buf, vs_buf, "\t ,;\n");
        j=0;
        _snp_name.push_back(vs_buf[++j]);
        _chr.push_back(atoi(vs_buf[++j].c_str()));
        _bp.push_back(atoi(vs_buf[++j].c_str()));
    }
    i_snp_info.close();
    _snp_num=_snp_name.size();
    cout<<_snp_num<<" SNPs specified in ["+snp_info_file+"]."<<endl;

    // save PLINK map file
    string map_file=_out+".map";
    ofstream omap(map_file.c_str());
    cout<<"Saving PLINK MAP file ..."<<endl;
    if(!omap) throw("Error: can not open the file ["+map_file+"] to write.");
    for(i=0; i<_snp_num; i++) omap<<_chr[i]<<"\t"<<_snp_name[i]<<"\t0\t"<<_bp[i]<<endl;
    omap.close();
    cout<<"PLINK MAP file has been saved in ["+map_file+"].\n"<<endl;

    // read the filenames of the raw genotype data files
    ifstream i_fnames(fname_file.c_str());
    if(!i_fnames) throw("Error: can not open the file ["+fname_file+"] to read.");
    vector<string> fnames;
    while(getline(i_fnames, str_buf)){
        if(StrFunc::split_string(str_buf, vs_buf)==1) fnames.push_back(vs_buf[0]);
    }
    i_fnames.close();
    _indi_num=fnames.size();
    cout<<_indi_num<<" raw genotype data filenames specified in ["+fname_file+"]."<<endl;

    // read raw genotype file
    _snp_1.resize(_snp_num);
    _snp_2.resize(_snp_num);
    for(i=0; i<_snp_num; i++){
        _snp_1[i].reserve(_indi_num);
        _snp_2[i].reserve(_indi_num);
    }
    cout<<"Reading the raw genotype files and saving the genotype data in PLINK PED format ..."<<endl;
    cout<<"(SNP genotypes with GenCall rate < "<<GC_cutoff<<" are regarded as missing)"<<endl;
    string ped_file=_out+".ped";
    ofstream oped(ped_file.c_str());
    if(!oped) throw("Error: can not open the file ["+ped_file+"] to read.");
    for(i=0; i<_indi_num; i++){
        read_one_IRG(oped, i, fnames[i], GC_cutoff);
        cout<<i+1<<" of "<<_indi_num<<" files.\r";
    }
    cout<<"Genotype data for "<<_indi_num<<" individuals have been save in the file ["+ped_file+"]."<<endl;
    oped.close();
}

char gcta::flip_allele(char a)
{
    if(a=='A') return('T');
    else if(a=='C') return('G');
    else if(a=='G') return('C');
    else if(a=='T') return('A');
    else return('0');
}

void gcta::read_one_IRG(ofstream &oped, int ind, string IRG_fname, double GC_cutoff)
{
    ifstream i_IRG(IRG_fname.c_str());
    if(!i_IRG) throw("Error: can not open the file ["+IRG_fname+"] to read.");
    char a1='0', a2='0';
    int i=0, j=0, snp_num=0;
    double GC=0.0;
    string str_buf, fid, pid;
    vector<string> vs_buf;
    for(i=0; i<2; i++) getline(i_IRG, str_buf);
    StrFunc::split_string(str_buf, vs_buf);
    bool oldversion=false;
    if(vs_buf[2]=="1.1.9") oldversion=true;
    for(i=0; i<3; i++) getline(i_IRG, str_buf);
    StrFunc::split_string(str_buf, vs_buf);
    snp_num=atoi(vs_buf[2].c_str());
    if(snp_num!=_snp_num) throw("Error: the number of SNPs specified in the summary file does not match that in the raw genotype file ["+IRG_fname+"].");
    for(i=0; i<6; i++) getline(i_IRG, str_buf);
    for(i=0; i<snp_num; i++){
        getline(i_IRG, str_buf);
        StrFunc::split_string(str_buf, vs_buf, "\t ,;\n");
        if(vs_buf[0]!=_snp_name[i]) throw("Error: the SNP ["+vs_buf[0]+"] specified in the summary file does not match that in the raw genotype file ["+IRG_fname+"]. The order of the SNPs have been changed?");
        pid=vs_buf[1];
        if(oldversion) fid=pid;
        else fid=vs_buf[2];
        if(oldversion){
            a1=vs_buf[2][0];
            a2=vs_buf[3][0];
            GC=atof(vs_buf[8].c_str());
        }
        else{
            GC=atof(vs_buf[3].c_str());
            a1=vs_buf[6][0];
            a2=vs_buf[7][0];
        }
        if(GC<GC_cutoff) a1=a2='0';
        if(i==0) oped<<fid<<" "<<pid<<" -9 -9 -9 -9 ";
        oped<<a1<<" "<<a2<<" ";
    }
    oped<<endl;
    i_IRG.close();
}

void gcta::save_plink()
{
    save_famfile();
    save_bimfile();
    save_bedfile();
}

void gcta::save_bedfile()
{
    int i=0, pos=0, j=0;
    string OutBedFile=_out+".bed";
	fstream OutBed(OutBedFile.c_str(), ios::out|ios::binary);
	if(!OutBed) throw("Error: can not open the file ["+OutBedFile+"] to write.");
	cout<<"Writing genotypes to PLINK BED file ["+OutBedFile+"] ..."<<endl;
	bitset<8> b;
	char ch[1];
    b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    b.reset();
    b.set(0);  b.set(1);  b.set(3);  b.set(4);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    b.reset();
    b.set(0);
    ch[0] = (char)b.to_ulong();
    OutBed.write(ch,1);
    for(i=0; i<_include.size(); i++){
        pos=0;
        b.reset();
        for(j=0; j<_keep.size(); j++){
            b[pos++]=(!_snp_2[i][j]);
            b[pos++]=(!_snp_1[i][j]);
            if(pos>7 || j==_keep.size()-1){
                ch[0]=(char)b.to_ulong();
                OutBed.write(ch,1);
                pos=0;
                b.reset();
            }
        }
    }
    OutBed.close();
}

void gcta::save_famfile()
{
    string famfile=_out+".fam";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the fam file "+famfile+" to save!");
	cout<<"Writing PLINK FAM file to ["+famfile+"] ..."<<endl;
	int i=0;
	for(i=0; i<_keep.size(); i++){
		Fam<<_fid[_keep[i]]<<"\t"<<_pid[_keep[i]]<<"\t"<<_fa_id[_keep[i]]<<"\t"<<_mo_id[_keep[i]]<<"\t"<<_sex[_keep[i]]<<"\t"<<_pheno[_keep[i]]<<endl;
	}
	Fam.close();
	cout<<_keep.size()<<" individuals to be saved to ["+famfile+"]."<<endl;
}

void gcta::save_bimfile()
{
	int i=0;
	string bimfile=_out+".bim";
	ofstream Bim(bimfile.c_str());
	if(!Bim) throw("Error: can not open the file ["+bimfile+"] to write.");
	cout<<"Writing PLINK bim file to ["+bimfile+"] ..."<<endl;
	for(i=0; i<_include.size(); i++){
		Bim<<_chr[_include[i]]<<"\t"<<_snp_name[_include[i]]<<"\t"<<_genet_dst[_include[i]]<<"\t"<<_bp[_include[i]]<<"\t"<<_allele1[_include[i]]<<"\t"<<_allele2[_include[i]]<<endl;
	}
	Bim.close();
	cout<<_include.size()<<" SNPs to be saved to ["+bimfile+"]."<<endl;
}
