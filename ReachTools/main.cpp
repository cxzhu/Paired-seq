//
// Chenxu Zhu (chz272@ucsd.edu)
// 6/29/2018
//

#include <iostream>
#include "reachtools.h"
#include "readgenome.h"
#include "cxstring.h"
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <unistd.h>

using namespace std;

void convert_reads(string read1, string read2);
void convert_reads_help();
void convRead2(string sample);
void convRead2_help();
void extractValidBC(string read);
void extractValidBC_help();
void splitBam(string bam);
void splitBam_help();
void bam2Mtx(string bam, string ref, string name, int cutoff);
void bam2Mtx2(string bam, string ref);
void bam2Mtx_help();
void rmdup(string bam);
void rmdup_help();
void help();
void pairBarcodes(string bam, string target);
void pairBarcodes_help();
void pairBC2R1(string bam, string target);
void pairBC2R1_help();
void pre_process(string read2);
void combine(string sample);
void pre_process_help();
void combine_help();
void convert(string sample);
void convert_help();
/*
void pca(string string);
void pca_help();

*/

string version = "2018.10.24";




int main(int argc, char** argv) {
	if(argc < 2){
		help();
		return 1;
	}
	string mod(argv[1]);
	if(mod == "-h" || mod == "help" || mod == "--help"){
		help();
		return 0;
	}

	//convert reads, convert R2 to cell barcode and UMI, merge R2 to readname of R1;
   if(mod == "convertReads"){
		if(argc < 4){
			convert_reads_help();
			return 1;
		}
		convert_reads(argv[2], argv[3]);
		return 0;
	}

	if(mod == "splitBam"){
		if(argc < 3){
			splitBam_help();
			return 1;
		}
		splitBam(argv[2]);
		return 0;
	}

	if(mod == "extVal"){
		if(argc < 3){
			extractValidBC_help();
			return 1;
		}
		extractValidBC(argv[2]);
		return 0;
	}

	if(mod == "bam2Mtx"){
		if(argc < 3){
			bam2Mtx_help();
			return 1;
		}
		int cutoff = 5;
		if(argc == 5){
			cutoff = cxstring::str2int(argv[4]);

		}
		string ref = "/projects/ps-renlab/chz272/annotations/human_hg19_genomic_region_annotation/hg19.gencode.v19.txt";
		// A reference file for DNA bins/peaks or RNA genes is needed here
		//DNA
		//chrom pos_start pos_end bin_name_1 bin_name_2
		//chr1 0 1000 chr1:0-1000 chr1:0-1000
		//chr1 1000 2000 chr1:1000-2000 chr1:1000-2000
		//
		// a representative mm10 in 5kb-bin ref file is in https://github.com/cxzhu/Paired-Tag/blob/master/refereces/mm10.bin5k.txt.zip

		//RNA
		//chrom pos_start pos_end gene_ID gene_symbol
		//chr1 3072253 3075322 ENSMUSG00000102693 4933401J01Rik
		//chr1 3101016 3103125 ENSMUSG00000064842 Gm26206
		// 
		// for mm10, you can access it from: https://github.com/cxzhu/Paired-Tag/blob/master/refereces/mm10.RNA.txt
		
		if(argc >= 4){
			ref = argv[3];
			if(ref == "mm10"){
				ref = "/projects/ps-renlab/chz272/annotations/mm10.big.txt";
			}
			else if(ref == "mm10m"){
				ref = "/projects/ps-renlab/chz272/annotations/mm10.big.backup";
			}
			else if(ref == "mm10utr"){
				ref = "/projects/ps-renlab/chz272/genome_ref/utr3_mm10.xls";
			}
			else if(ref == "hg19bin"){
				ref = "/projects/ps-renlab/chz272/genome_ref/hg19.bin10k.txt";
			}
			else if(ref == "mm10_10k"){
				ref = "/projects/ps-renlab/chz272/genome_ref/mm10.bin10k.txt";
			}
			else if(ref == "mm10_1k"){
				ref = "/projects/ps-renlab/chz272/genome_ref/mm10.bin1k.txt";
			}
			else if(ref == "mm10_reg"){
				ref= "/projects/ps-renlab/chz272/genome_ref/refdata-cellranger-atac-mm10-1.0.1/genes/Regulatory.xls";
			}
			else if(ref == "hg19"){
				ref="/projects/ps-renlab/chz272/annotations/human_hg19_genomic_region_annotation/hg19.gcode.v19.txt";
			}
			else if(ref == "hg19_1k"){
				ref="/projects/ps-renlab/chz272/genome_ref/hg19.bin1k.txt";
			}

		}
		bam2Mtx(argv[2], ref, argv[3], cutoff) ;
		return 0;
	}



	if(mod == "convRead2"){
		if(argc < 3){
			convRead2_help();
			return 1;
		}
		convRead2(argv[2]);
		return 0;
	}

	if(mod == "rmdup"){
		if(argc < 3){
			rmdup_help();
			return 1;
		}
		rmdup(argv[2]);
		return 0;
	}

	if(mod == "pre_process"){
		if(argc < 3){
			pre_process_help();
			return 1;
		}
		pre_process(argv[2]);
		return 0;
	}

	if(mod == "combine"){
		if(argc < 3){
			combine_help();
			return 1;
		}
		combine(argv[2]);
		return 0;
	}

	if(mod == "convert"){
		if(argc < 3){
			convert_help();
			return 1;
		}
		convert(argv[2]);
		return 0;
	}

	if(mod == "pairBarcodes"){
		if(argc < 4){
			pairBarcodes_help();
			return 1;
		}
		pairBarcodes(argv[2], argv[3]);
		return 0;
	}
	if(mod == "pairBC2R1"){
		if(argc < 4){
			pairBarcodes_help();
			return 1;
		}
		pairBarcodes(argv[2], argv[3]);
		return 0;
	}

	return 0;
}

//  local functions
void help(){
	cout << "ReachTools, processiong of REACH/REChIP-seq data." << "\nVersion: " << version << " (cxzhu@pku.edu.cn)" << endl;
	cout << endl << "\tconvertReads\tExtract cell barcodes from read2, filt low quality reads and merge to 1 raw fastq file." << endl;
	cout << "\tconvRead2\tExtract cell barcodes from read2, filt low quality reads and merge to 1 raw fastq file. \n\t\t\tNew version for tandom N in 5' of read2." << endl;
	cout << "\tpre_processs\tExtract UMI and BC sequence from read2." << endl;
	cout << "\tpairBarcodes\tPair the BC1-3(from pre_proc), BC4 and UMI to target sam/bam file." << endl;
//	cout << "\textVal\t\tExtract reads with valid barcodes from merged read file." << endl;
	cout << "\tsplitBam\tSpit bam file into FILE_DNA.bam, FILE_RNA.bam and FILE_UND.bam." << endl;
	cout << "\trmdup\t\tRemove duplicated reads based on UMI and mapped position." << endl;
	cout << "\tbam2Mtx\t\tTransform splitted bam file into Matrix File." << endl;
//	cout << endl << "\tgenerateMatrix\t\tGenerate data-matrix from peak calling file and bam file." << endl;

	return;
}

void convert_reads_help(){
	cout << "reachtools convertReads R1.fq.gz R2.fq.gz" << endl;
	return;
}

void convert_reads(string read1, string read2){


	/// file type
	vector<string> tmp;
	tmp = cxstring::split(read1, ".");
	string s1 = "cat ";
	if(tmp[tmp.size() - 1] == "gz")s1 = "zcat ";
	if(tmp[tmp.size() - 1] == "bz2")s1 = "bzcat ";

	{
		///if test
		ifstream iftest;
		iftest.open(read1);
		if(!iftest){
			cout << "Cannot open read1 file\n";
			exit(-1);
		}
		iftest.close();
		iftest.open(read2);
		if(!iftest){
			cout << "Cannot open read2 file\n";
			exit(-1);
		}
		iftest.close();
	}

	string s2 = read1;
	string s3 = s1 + s2;


	FILE * in_read1;
	in_read1 = popen(s3.c_str(), "r");

	s2 = read2;
	s3 = s1 + s2;

	FILE * in_read2;
	in_read2 = popen(s3.c_str(), "r");

	tmp = cxstring::split(read1, "_");
	s3 = "gzip - > " + tmp[0] + "_BC.fq.gz";
	FILE *out_file;
	out_file = popen(s3.c_str(), "w");

	char buffer[1000];
	int line_number = 0;
	int _line_number = 0;

	//read fastq file
	line_number = 0;
	int mega_line_number = 0;

	int num_length = 0;
	int num_docking = 0;
	int num_bcQual = 0;
	int num_bcValid = 0;
	int num_dna = 0;
	int num_rna = 0;
	while(fgets(buffer, sizeof(buffer), in_read1)) {
		//process each lines

		++line_number;
		++_line_number;
		if (line_number == 10000000) {
			++mega_line_number;
			line_number = 0;
			cout << mega_line_number << "0 M reads processed..." << endl;
		}

		// initial vars for pair-end reads
		string r1_read_name;
		string r1_sequence;
		string r1_mark;
		string r1_quality;
		string r2_read_name;
		string r2_sequence;
		string r2_mark;
		string r2_quality;

//	  stringstream tmp(buffer);
		stringstream tmp;
		tmp << buffer;
		tmp >> r1_read_name;
		r1_read_name = cxstring::chomp(r1_read_name);
		// read read1
		tmp.str("");
		fgets(buffer, sizeof(buffer), in_read1);
		tmp << buffer;
		tmp >> r1_sequence;
		fgets(buffer, sizeof(buffer), in_read1);
		tmp << buffer;
		tmp >> r1_mark;
		fgets(buffer, sizeof(buffer), in_read1);
		tmp << buffer;
		tmp >> r1_quality;

		// read read2
		tmp.str("");
		fgets(buffer, sizeof(buffer), in_read2);
		tmp << buffer;
		tmp >> r2_read_name;
		tmp.str("");
		fgets(buffer, sizeof(buffer), in_read2);
		tmp << buffer;
		tmp >> r2_sequence;
		fgets(buffer, sizeof(buffer), in_read2);
		tmp << buffer;
		tmp >> r2_mark;
		fgets(buffer, sizeof(buffer), in_read2);
		tmp << buffer;
		tmp >> r2_quality;
		tmp.str("");

		//transfer R2 to barcode infor
		string read_id;
		if(r2_sequence.length() < 128)continue;
		++num_length;

		//docking barcodes from Read2
		reachtools::docking(r2_sequence, r2_quality);
		if(r2_sequence == "Low Quality")continue;
		++num_docking;
		if(r2_quality == "Low Quality")continue;
		++num_bcQual;

		//extract barcodes
		read_id = reachtools::extract_barcode(r2_sequence); //, bc_lib);
		vector<string> tm;
		tm = cxstring::split(read_id, ":");
		if(tm[0] == "d"){
			++num_dna;
		}
		else if(tm[0] == "r"){
			++num_rna;
		}
		if(tm[1] != "00" && tm[2] != "00" && tm[3] != "00" && tm[4] != "00")++num_bcValid;
		r1_read_name = r1_read_name + ':' + read_id;
		fputs((r1_read_name + "\n").c_str(), out_file);
		fputs((r1_sequence + "\n").c_str(), out_file);
		fputs((r1_mark + "\n").c_str(), out_file);
		fputs((r1_quality + "\n").c_str(), out_file);
	}
	pclose(in_read1);
	pclose(in_read2);
	pclose(out_file);

	s3 = tmp[0] + "_BC_extraction_report.txt";
	ofstream of;
	of.open(s3);
	of << _line_number << " total reads processed." << endl;
	of << "\t" << num_length << " reads have length enough." << endl;
	of << "\t" << num_docking << " reads have valid docking sequence." << endl;
	of << "\t" << num_bcQual << " reads have all barcodes pass Q20." << endl;
	of << "\t" << num_bcValid << " reads have all barcodes valid." << endl;
	of << "\t" << num_dna << " DNA reads and " << num_rna << " RNA reads." << endl;
	of << "Reads with non-valid barcodes were reserved." << endl;
	of.close();
	return;
}

void extractValidBC_help(){
	cout << "reachtools extVal sample_BC.fq.gz" << endl;
	return;

}
void extractValidBC(string read){

	return;
}


void splitBam_help(){
	cout << "reachtools splitBam in.bam/in.sam" << endl;
	return;

}
void splitBam(string bam){
	string s1 = "samtools view -h ";
	{
		///if test
		ifstream iftest;
		iftest.open(bam);
		if(!iftest){
			cout << "Cannot open read1 file\n";
			exit(-1);
		}
		iftest.close();
	}

	string s2 = bam;
	string s3 = s1 + s2;

	FILE * in_file;
	in_file = popen(s3.c_str(), "r");
	FILE * dn_file;
	s1 = "cat > ";
	s2 = bam.substr(0, bam.length()-4) + "_DNA.tmp";
	s3 = s1 + s2;
	dn_file = popen(s3.c_str(), "w");
	FILE * rn_file;
	s1 = "cat > ";
	s2 = bam.substr(0, bam.length()-4) + "_RNA.tmp";
	s3 = s1 + s2;
	rn_file = popen(s3.c_str(), "w");
	FILE * nn_file;
	s1 = "cat > ";
	s2 = bam.substr(0, bam.length()-4) + "_UND.tmp";
	s3 = s1 + s2;
	nn_file = popen(s3.c_str(), "w");

	char buffer[2000];
	while(fgets(buffer, sizeof(buffer), in_file)){
		string line(buffer);
		cxstring::chomp(line);
		//header
		if(line.substr(0, 1) == "@"){
			fputs(line.c_str(), dn_file);
			fputs(line.c_str(), rn_file);
			fputs(line.c_str(), nn_file);
		}
		else{
			vector<string> tmp;
			tmp = cxstring::split(line, ":");
			if(tmp[7] == "d"){
				fputs(line.c_str(), dn_file);
			}
			else if(tmp[7] == "r"){
				fputs(line.c_str(), rn_file);
			}
			else{
				fputs(line.c_str(), nn_file);
			}
		}
	}
	pclose(in_file);
	pclose(dn_file);
	pclose(rn_file);
	pclose(nn_file);
	s1 = "samtools view -bS ";
	s2 = bam.substr(0, bam.length()-4) + "_DNA.tmp";
	s3 = s1 + s2 + ">";
	s2 = bam.substr(0, bam.length()-4) + "_DNA.bam";
	s3 = s3 + s2;
	system(s3.c_str());
	s2 = bam.substr(0, bam.length()-4) + "_DNA.tmp";
	s3 = "rm " + s2;
	system(s3.c_str());
	s1 = "samtools view -bS ";
	s2 = bam.substr(0, bam.length()-4) + "_RNA.tmp";
	s3 = s1 + s2 + ">";
	s2 = bam.substr(0, bam.length()-4) + "_RNA.bam";
	s3 = s3 + s2;
	system(s3.c_str());
	s2 = bam.substr(0, bam.length()-4) + "_RNA.tmp";
	s3 = "rm " + s2;
	system(s3.c_str());
	s1 = "samtools view -bS ";
	s2 = bam.substr(0, bam.length()-4) + "_UND.tmp";
	s3 = s1 + s2 + ">";
	s2 = bam.substr(0, bam.length()-4) + "_UND.bam";
	s3 = s3 + s2;
	system(s3.c_str());
	s2 = bam.substr(0, bam.length()-4) + "_UND.tmp";
	s3 = "rm " + s2;
	system(s3.c_str());
	return;
}


void bam2Mtx_help(){
	cout << "reachtools Bam2Mtx in.bam/in.sam [ref]\n\tref is optional for RNA-seq but mandatory for ATAC-seq" << endl;
	return;
}

void bam2Mtx(string bam, string ref, string nam, int cutoff){
	// ref ref annotation
	string s1 = "cat ";
	string s2 = ref;
	string s3 = s1 + s2;
	{
		//if test
				ifstream iftest;
				iftest.open(ref);
				if(!iftest){
					cout << "Cannot open ref file\n";
					exit(-1);
				}
				iftest.close();
	}
	FILE * in_file;
	in_file = popen(s3.c_str(), "r");
	char buffer[2000];
	cout << "Reading ref annotation..." << endl;
	map<string, map <int, string>> genes;
	map<string, map <int, string>>::iterator igenes;
	map<int, string>::iterator iigenes;
	while(fgets(buffer, sizeof(buffer), in_file)){

		string line(buffer);
		line = cxstring::chomp(line);
		vector<string> tmp = cxstring::split(line, "\t");
		/*
		if(tmp.size() < 4){
			cout << tmp.size() << endl;
			cout << line << endl;
			exit(0);
		}
		*/
		string chr = tmp[0];
		int pss = cxstring::str2int(tmp[1]);
		int pse = cxstring::str2int(tmp[2]);
		string gene_id = tmp[3];
		for(int pos = int(pss/1000); pos <= int(pse/1000); pos++){
			string gene_line = tmp[1] + "\t" + tmp[2] + "\t" + gene_id;
			genes[chr][pos] = gene_line;
		}
	} 
	pclose(in_file);


	//read sam file
	s1 = "samtools view ";
	{
		///if test
		ifstream iftest;
		iftest.open(bam);
		if(!iftest){
			cout << "Cannot open read1 file\n";
			exit(-1);
		}
		iftest.close();
	}

	s2 = bam;
	s3 = s1 + s2;

	in_file = popen(s3.c_str(), "r");
	//int count = 0;
	//int mega_count = 0;
	cout << "Processing sam/bam file..." << endl;
	// map[gene][cell_id][umi] hash
	// map[gene] genelist
	// map[cell_id] celllist
	map<string, map<string, map<string, bool>>> hash;
	map<string, map<string, map<string, bool>>>::iterator ihash;
	map<string, map<string, bool>>::iterator iihash;
	map<string, bool>::iterator iiihash;
	map<string, int> genelist;
	//map<string, int>::iterator igenelist;
	map<string, int> celllist;
	//map<string, int>::iterator icelllist;
	while(fgets(buffer, sizeof(buffer), in_file)){
		/*
		++count;
		if(count > 10000000){
			count = 0;
			++mega_count;
			cout << "\t" << mega_count << "0M reads read..." << endl;
		}
		*/
		string line(buffer);
		samline sl;
		sl.init(line);
		// check chr exists
		igenes=genes.find(sl.chr);
		if(igenes==genes.end()){
			continue;
		}
		int spos = int(sl.pos/1000);
		iigenes=igenes->second.find(spos);
		if(iigenes==igenes->second.end()){
			continue;
		}
		string gene_line = iigenes->second;
		vector<string> tmp = cxstring::split(gene_line, "\t");
		int pss = cxstring::str2int(tmp[0]);
		if(sl.pos < pss){
			continue;
		}
		int pse = cxstring::str2int(tmp[1]);
		if(sl.pos > pse){
			continue;
		}
		string gene_id = tmp[2];
		/*
		vector<string> sp;
		sp = cxstring::split(sl.readname, ":");
		int l = sp.size();
		string cell_id = sp[l-5] + ":" + sp[l-4] + ":" + sp[l-3] + ":" + sp[l-2];
		string umi = sp[l-1];
		50:77:25:08:CCTTTATATA
		*/
		string cell_id = sl.readname.substr(sl.readname.length()-22, 11);
		string umi = sl.readname.substr(sl.readname.length()-10, 10);
		//hash[gene_id][cell_id][umi] = true;
		hash[cell_id][gene_id][umi] = true;
		genelist[gene_id] = 0;
		celllist[cell_id] = 0;
	}
	pclose(in_file);
	cout << "Passing cell post filter to new hash..." << endl;
	// cutoff = 50 reads per cell
	//
	//
	//
	////
	//
	int total = 0;
	map<string, map<string, int>> new_hash;
	map<string, map<string, int>>::iterator inew_hash;
	map<string, int>::iterator iinew_hash;
  map<string, int> new_genelist;
  map<string, int>::iterator inew_genelist;
  map<string, int> new_celllist;
  map<string, int>::iterator inew_celllist;
	map<string, bool> empty_iimap;
	map<string, map<string, bool>> empty_imap;
	map<string, map<string, map<string, bool>>> empty_map;
	for(ihash = hash.begin(); ihash != hash.end(); ihash++){ // ihash->first = cell_id, ihash->second = map<genelist, umi>
		//per cell
		int n_cell_umi = 0;
		map<string, int> cur_genelist;
		map<string, int>::iterator i_cur_genelist;
		for(iihash = ihash->second.begin(); iihash != ihash->second.end(); iihash++){ // iihash->first = genelist, iihash->second = <umi, T>
			//per gene
			int n_umi = iihash->second.size(); 
			//int n_umi = 0;
			/*	
			for(iiihash = iihash->second.begin(); iiihash != iihash->second.end(); iiihash++){ 
				++n_umi;
				//iihash->second.swap(empty_iimap);
				//iihash->second.clear();
			}
			*/	
			cur_genelist[iihash->first] = n_umi;
			n_cell_umi += n_umi;
			iihash->second.swap(empty_iimap);
			iihash->second.clear();
		}
		if(n_cell_umi >= cutoff){ ////// SETUP the CUTOFF!!
			for(i_cur_genelist = cur_genelist.begin(); i_cur_genelist != cur_genelist.end(); ++i_cur_genelist){
				++total;
				new_hash[i_cur_genelist->first][ihash->first] = i_cur_genelist->second;
				new_genelist[i_cur_genelist->first] = 1;
			}
			new_celllist[ihash->first] = 1;
			//cout << total << endl; ///////////////////////////////////////////////// DIAG
		}
		ihash->second.swap(empty_imap);
		ihash->second.clear();
	}
	// release memory
	cout << "Cleaning memory..." << endl;
	hash.clear();
	celllist.clear();
	genelist.clear();

		
	// output new mtx file
	cout << "Outputing MatrixFactory Format File..." << endl;
	/*
	FILE * dn_file;
	s1 = "cat > ";
	s2 = bam.substr(0, bam.length()-4) + "_DNA.tmp";
	s3 = s1 + s2;
	dn_file = popen(s3.c_str(), "w");
	fputs(line.c_str(), dn_file);
	*/

	string name = bam.substr(0, 12);
	name = name + "_" + nam;
	s1 = "mkdir ";
	s2 = name + "_mtx";
	s3 = s1 + s2;
	if(access(s2.c_str(), 0) == -1){
		system(s3.c_str());
	}

	s1 = "cat > ";
	s2 = name + "_mtx/genes.tsv";
	s3 = s1 + s2;
	FILE *  out_file;
	out_file = popen(s3.c_str(), "w");
	int gene_order = 0;
	for(inew_genelist = new_genelist.begin(); inew_genelist != new_genelist.end(); inew_genelist++){
		++gene_order;
		vector<string> tmp;
		string tstring = inew_genelist->first;
		tmp = cxstring::split(tstring, " ");
		string output = tmp[0] + "\t" + tmp[1] + "\n";
		//cout << output; ////////////////////////////////////////////////////////////////////////// DIAG
		fputs(output.c_str(), out_file);
		inew_genelist->second = gene_order;
	} 
	pclose(out_file);

	s1 = "cat > ";
	s2 = name + "_mtx/barcodes.tsv";
	s3 = s1 + s2;
	out_file = popen(s3.c_str(), "w");
	int cell_order = 0;
	for(inew_celllist = new_celllist.begin(); inew_celllist != new_celllist.end(); inew_celllist++){
		++cell_order;
		string output = inew_celllist->first + "\n";
		fputs(output.c_str(), out_file);
		inew_celllist->second = cell_order;
	}
	pclose(out_file);

	s1 = "cat > ";
	s2 = name + "_mtx/matrix.mtx";
	s3 = s1 + s2;
	out_file = popen(s3.c_str(), "w");
	fputs("%%MatrixMarket matrix coordinate real general\n", out_file);
	fputs("%\n", out_file);
	string output;
	output = cxstring::int2str(gene_order) + " " + cxstring::int2str(cell_order) + " " + cxstring::int2str(total) + "\n";
	//output = cxstring::int2str(new_genelist.size()) + " " + cxstring::int2str(new_celllist.size()) + " " + cxstring::int2str(total) + "\n";
	fputs(output.c_str(), out_file);

	for(inew_hash = new_hash.begin(); inew_hash != new_hash.end(); inew_hash++){
		for(iinew_hash = inew_hash->second.begin(); iinew_hash != inew_hash->second.end(); iinew_hash++){
			//string output = inew_hash->first + " " + iinew_hash->first + " " + cxstring::int2str(iinew_hash->second) + "\n";
			string output;
			inew_genelist = new_genelist.find(inew_hash->first);
			inew_celllist = new_celllist.find(iinew_hash->first);
			output = cxstring::int2str(inew_genelist->second) + " " + cxstring::int2str(inew_celllist->second) + " " + cxstring::int2str(iinew_hash->second) + "\n";
			fputs(output.c_str(), out_file);
		}
	}
	pclose(out_file);
	return;
}

void bam2Mtx2(string bam, string ref){
	// ref ref annotation
	string s1 = "cat ";
	string s2 = ref;
	string s3 = s1 + s2;
	{
		//if test
				ifstream iftest;
				iftest.open(ref);
				if(!iftest){
					cout << "Cannot open ref file\n";
					exit(-1);
				}
				iftest.close();
	}
	FILE * in_file;
	in_file = popen(s3.c_str(), "r");
	char buffer[2000];
	cout << "Reading ref annotation..." << endl;
	
	/* version 2018, using displaceable annotation
	map<string, map <int, string>> genes;
	map<string, map <int, string>>::iterator igenes;
	map<int, string>::iterator iigenes;
	while(fgets(buffer, sizeof(buffer), in_file)){

		string line(buffer);
		line = cxstring::chomp(line);
		vector<string> tmp = cxstring::split(line, "\t");

		string chr = tmp[0];
		int pss = cxstring::str2int(tmp[1]);
		int pse = cxstring::str2int(tmp[2]);
		string gene_id = tmp[3];
		for(int pos = int(pss/1000); pos <= int(pse/1000); pos++){
			string gene_line = tmp[1] + "\t" + tmp[2] + "\t" + gene_id;
			genes[chr][pos] = gene_line;
		}
	}
	*/
	// version 2019, try to update with un displaceable annotation
	map<string, map <int, map<string, int>>> genes;
	map<string, map <int, map<string, int>>>::iterator igenes;
	map<int, map<string, int>>::iterator iigenes;
	map<string, int>::iterator iiigenes;
	while(fgets(buffer, sizeof(buffer), in_file)){
		string line(buffer);
		line = cxstring::chomp(line);
		vector<string> tmp = cxstring::split(line, "\t");

		string chr = tmp[0];
		int pss = cxstring::str2int(tmp[1]);
		int pse = cxstring::str2int(tmp[2]);
		string gene_id = tmp[3];
		for(int pos = int(pss/1000); pos <= int(pse/1000); pos++){
			string gene_line = tmp[1] + "\t" + tmp[2] + "\t" + gene_id;
			genes[chr][pos][gene_line] = 1;
		}
	}
	pclose(in_file);
	


	//read sam file
	s1 = "samtools view ";
	{
		///if test
		ifstream iftest;
		iftest.open(bam);
		if(!iftest){
			cout << "Cannot open read1 file\n";
			exit(-1);
		}
		iftest.close();
	}

	s2 = bam;
	s3 = s1 + s2;

	in_file = popen(s3.c_str(), "r");
	//int count = 0;
	//int mega_count = 0;
	cout << "Processing sam/bam file..." << endl;
	// map[gene][cell_id][umi] hash
	// map[gene] genelist
	// map[cell_id] celllist
	map<string, map<string, map<string, bool>>> hash;
	map<string, map<string, map<string, bool>>>::iterator ihash;
	map<string, map<string, bool>>::iterator iihash;
	map<string, bool>::iterator iiihash;
	map<string, int> genelist;
	//map<string, int>::iterator igenelist;
	map<string, int> celllist;
	//map<string, int>::iterator icelllist;
	while(fgets(buffer, sizeof(buffer), in_file)){
		/*
		++count;
		if(count > 10000000){
			count = 0;
			++mega_count;
			cout << "\t" << mega_count << "0M reads read..." << endl;
		}
		*/
		string line(buffer);
		samline sl;
		sl.init(line);
		// check chr exists
		igenes=genes.find(sl.chr);
		if(igenes==genes.end()){
			continue;
		}
		int spos = int(sl.pos/1000);
		iigenes=igenes->second.find(spos);
		if(iigenes==igenes->second.end()){
			continue;
		}
		/*
		string gene_line = iigenes->second;
		/// insert new code here
		vector<string> tmp = cxstring::split(gene_line, "\t");
		int pss = cxstring::str2int(tmp[0]);
		if(sl.pos < pss){
			continue;
		}
		int pse = cxstring::str2int(tmp[1]);
		if(sl.pos > pse){
			continue;
		}
		string gene_id = tmp[2];
		string cell_id = sl.readname.substr(sl.readname.length()-22, 8);
		string umi = sl.readname.substr(sl.readname.length()-10, 10);
		//hash[gene_id][cell_id][umi] = true;
		hash[cell_id][gene_id][umi] = true;
		genelist[gene_id] = 0;
		celllist[cell_id] = 0;
		*/
		for(iiigenes=iigenes->second.begin(); iiigenes!=iigenes->second.end(); ++iiigenes){
			string gene_line = iiigenes->first;
			vector<string> tmp = cxstring::split(gene_line, "\t");
			int pss = cxstring::str2int(tmp[0]);
			if(sl.pos < pss){
				continue;
			}
			int pse = cxstring::str2int(tmp[1]);
			if(sl.pos > pse){
				continue;
			}
			string gene_id = tmp[2];
			string cell_id = sl.readname.substr(sl.readname.length()-22, 8);
			string umi = sl.readname.substr(sl.readname.length()-10, 10);
			//hash[gene_id][cell_id][umi] = true;
			hash[cell_id][gene_id][umi] = true;
			genelist[gene_id] = 0;
			celllist[cell_id] = 0;
		}
	}
	pclose(in_file);
	cout << "Passing cell post filter to new hash..." << endl;
	// cutoff = 50 reads per cell
	//
	//
	//
	int cutoff = 5; ////
	//
	int total = 0;
	map<string, map<string, int>> new_hash;
	map<string, map<string, int>>::iterator inew_hash;
	map<string, int>::iterator iinew_hash;
  map<string, int> new_genelist;
  map<string, int>::iterator inew_genelist;
  map<string, int> new_celllist;
  map<string, int>::iterator inew_celllist;
	map<string, bool> empty_iimap;
	map<string, map<string, bool>> empty_imap;
	map<string, map<string, map<string, bool>>> empty_map;
	for(ihash = hash.begin(); ihash != hash.end(); ihash++){ // ihash->first = cell_id, ihash->second = map<genelist, umi>
		//per cell
		int n_cell_umi = 0;
		map<string, int> cur_genelist;
		map<string, int>::iterator i_cur_genelist;
		for(iihash = ihash->second.begin(); iihash != ihash->second.end(); iihash++){ // iihash->first = genelist, iihash->second = <umi, T>
			//per gene
			int n_umi = iihash->second.size(); 
			//int n_umi = 0;
			/*	
			for(iiihash = iihash->second.begin(); iiihash != iihash->second.end(); iiihash++){ 
				++n_umi;
				//iihash->second.swap(empty_iimap);
				//iihash->second.clear();
			}
			*/	
			cur_genelist[iihash->first] = n_umi;
			n_cell_umi += n_umi;
			iihash->second.swap(empty_iimap);
			iihash->second.clear();
		}
		if(n_cell_umi >= cutoff){ ////// SETUP the CUTOFF!!
			for(i_cur_genelist = cur_genelist.begin(); i_cur_genelist != cur_genelist.end(); ++i_cur_genelist){
				++total;
				new_hash[i_cur_genelist->first][ihash->first] = i_cur_genelist->second;
				new_genelist[i_cur_genelist->first] = 1;
			}
			new_celllist[ihash->first] = 1;
			//cout << total << endl; ///////////////////////////////////////////////// DIAG
		}
		ihash->second.swap(empty_imap);
		ihash->second.clear();
	}
	// release memory
	cout << "Cleaning memory..." << endl;
	hash.clear();
	celllist.clear();
	genelist.clear();

		
	// output new mtx file
	cout << "Outputing MatrixFactory Format File..." << endl;
	/*
	FILE * dn_file;
	s1 = "cat > ";
	s2 = bam.substr(0, bam.length()-4) + "_DNA.tmp";
	s3 = s1 + s2;
	dn_file = popen(s3.c_str(), "w");
	fputs(line.c_str(), dn_file);
	*/

	string name = bam.substr(0, bam.length()-4);
	s1 = "mkdir ";
	s2 = name + "_mtx2";
	s3 = s1 + s2;
	if(access(s2.c_str(), 0) == -1){
		system(s3.c_str());
	}

	s1 = "cat > ";
	s2 = name + "_mtx2/genes.tsv";
	s3 = s1 + s2;
	FILE *  out_file;
	out_file = popen(s3.c_str(), "w");
	int gene_order = 0;
	for(inew_genelist = new_genelist.begin(); inew_genelist != new_genelist.end(); inew_genelist++){
		++gene_order;
		vector<string> tmp;
		string tstring = inew_genelist->first;
		tmp = cxstring::split(tstring, " ");
		string output = tmp[0] + "\t" + tmp[1] + "\n";
		//cout << output; ////////////////////////////////////////////////////////////////////////// DIAG
		fputs(output.c_str(), out_file);
		inew_genelist->second = gene_order;
	} 
	pclose(out_file);

	s1 = "cat > ";
	s2 = name + "_mtx2/barcodes.tsv";
	s3 = s1 + s2;
	out_file = popen(s3.c_str(), "w");
	int cell_order = 0;
	for(inew_celllist = new_celllist.begin(); inew_celllist != new_celllist.end(); inew_celllist++){
		++cell_order;
		string output = inew_celllist->first + "\n";
		fputs(output.c_str(), out_file);
		inew_celllist->second = cell_order;
	}
	pclose(out_file);

	s1 = "cat > ";
	s2 = name + "_mtx2/matrix.mtx";
	s3 = s1 + s2;
	out_file = popen(s3.c_str(), "w");
	fputs("%%MatrixMarket matrix coordinate real general\n", out_file);
	fputs("%\n", out_file);
	string output;
	output = cxstring::int2str(gene_order) + " " + cxstring::int2str(cell_order) + " " + cxstring::int2str(total) + "\n";
	//output = cxstring::int2str(new_genelist.size()) + " " + cxstring::int2str(new_celllist.size()) + " " + cxstring::int2str(total) + "\n";
	fputs(output.c_str(), out_file);

	for(inew_hash = new_hash.begin(); inew_hash != new_hash.end(); inew_hash++){
		for(iinew_hash = inew_hash->second.begin(); iinew_hash != inew_hash->second.end(); iinew_hash++){
			//string output = inew_hash->first + " " + iinew_hash->first + " " + cxstring::int2str(iinew_hash->second) + "\n";
			string output;
			inew_genelist = new_genelist.find(inew_hash->first);
			inew_celllist = new_celllist.find(iinew_hash->first);
			output = cxstring::int2str(inew_genelist->second) + " " + cxstring::int2str(inew_celllist->second) + " " + cxstring::int2str(iinew_hash->second) + "\n";
			fputs(output.c_str(), out_file);
		}
	}
	pclose(out_file);
	return;
}


void convRead2(string sample){
	string s1 = "zcat ";
	string s2 = sample + "_R1.fq.gz";
	string s3 = s1 + s2;
	FILE * r1;
	r1 = popen(s3.c_str(), "r");
	FILE * r2;
	s2 = sample + "_R2.fq.gz";
	s3 = s1 + s2;
	r2 = popen(s3.c_str(), "r");
	FILE * outfile;
	s1 = "gzip - > ";
	s2 = sample + "_BC.fq.gz";
	s3 = s1 + s2;
	outfile = popen(s3.c_str(), "w");
	char buffer[1000];
	stringstream tmp;
	int tot = 0;
	int doc = 0;
	int val = 0;
	int dna = 0;
	int rna = 0;
	while(fgets(buffer, sizeof(buffer), r1)){
		string readname;
		tmp << buffer;
		tmp >> readname;
		tmp.str("");


		fqline r1_line;
		r1_line.read_part_record(r1, readname);
		
		fqline r2_line;
		r2_line.read_full_record(r2);

		bc_library lib;
		lib.init();
		read2 read_2;
		read_2.init(r2_line.seq);
		read_2.trim();
		//cout << read_2.dock << endl; //////////////////////////////////////////////////
		++tot;

		if(read_2.dock<4){
			//cerr << read_2.dock;
			continue;
		}
		++doc;


		read_2.extract_barcode(lib, 8, 0);
		if(read_2.type == "d")++dna;
		if(read_2.type == "r")++rna;

		if(read_2.is_valid())val++;
		r1_line.readname = r1_line.readname + ":" + read_2.type + ":" + read_2.bc + ":" + read_2.umi;

		
		r1_line.write_record(outfile);
	}
	cout << "BC extraction finished for " << sample << endl;
	cout << "Total " << tot << " reads" << endl;
	cout << "Docking " << doc << " reads" << endl;
	cout << "Valid " << val << " reads" << endl;
	cout << "DNA " << dna << " reads" << endl;
	cout << "RNA " << rna << " reads" << endl;
}
void convRead2_help(){
	cout << "reachtools convRead2 [sample_name]\n\tReads file should be sample_name_R1.fq.gz,sample_name_R2.fq.gz" << endl;
	return;
}

void rmdup(string bam){
	string name = bam.substr(0, bam.length()-4);
	string s1 = "samtools view ";
	string s2 = bam;
	string s3 = s1 + s2;
	FILE * inbam;
	inbam = popen(s3.c_str(), "r");
	map<string, map<int, map<string, bool>>> hash;
	char buffer[2000];
	while(fgets(buffer, sizeof(buffer), inbam)){
		string line(buffer);
		vector<string> tmp = cxstring::split(line, "\t");
		string chr = tmp[2];
		int pos = cxstring::str2int(tmp[3]);
		// 95:78:49:08:CCCCTCGTGG
		string umi = tmp[0].substr(tmp[0].length()-22, 22);
		hash[chr][pos][umi] = true;
	}
	pclose(inbam);

	s1 = "samtools view -h ";
	s3 = s1 + s2;
	inbam = popen(s3.c_str(), "r");
	s1 = "samtools view -b - > ";
	s2 = name + "_rmdup.bam";
	s3 = s1 + s2;
	FILE * outbam;
	outbam = popen(s3.c_str(), "w");

	while(fgets(buffer, sizeof(buffer), inbam)){
		string line(buffer);
		if(line.substr(0, 1) == "@"){
			fputs(line.c_str(), outbam);
		}	
		else{
			vector<string> tmp = cxstring::split(line, "\t");
			string chr = tmp[2];
			int pos = cxstring::str2int(tmp[3]);
			string umi = tmp[0].substr(tmp[0].length()-22, 22);
			
			if(hash[chr][pos][umi] == false)continue;

			fputs(line.c_str(), outbam);
			hash[chr][pos][umi] = false;
		}
	}
	pclose(inbam);
	pclose(outbam);

}
void rmdup_help(){
	cout << "reachtools rmdup input.bam\nOutput is input_rmdup.bam" << endl;
}


void pre_process(string r2){
	string s1 = "zcat ";
	string s2 = r2;
	string s3 = s1 + s2;
	FILE * infile;
	infile = popen(s3.c_str(), "r");
	s1 = "gzip - > ";
	s2 = r2.substr(0, r2.length()-6) + "_UMIBC1.fq.gz";
	s3 = s1 + s2;
	FILE * outfile;
	outfile = popen(s3.c_str(), "w");
	char buffer[2000];
	fqline in_line;
	while(fgets(buffer, sizeof(buffer), infile)){
		string line(buffer);
		line = cxstring::chomp(line);
		in_line.read_part_record(infile, line);
		read2 read_2;

		read_2.init(in_line.seq);
		read_2.trim();
		string new_seq = read_2.sbc1 + read_2.sbc2 + read_2.sbc3 + read_2.sbc4;
		string umi = read_2.umi;
		//string umi = in_line.seq.substr(0, 10);
		in_line.seq = new_seq;
		in_line.qual = in_line.qual.substr(0, in_line.seq.length());
		if(in_line.seq.length()!=24)continue;
		in_line.seq = new_seq;
		string a;
		stringstream as;
		as << in_line.readname;
		as >> a;
		in_line.readname = a + ":" + umi;
		//in_line.readname = in_line.readname + ":" + umi;
		in_line.write_record(outfile);
	}
	pclose(infile);
	pclose(outfile);

}

void combine(string r2){
	int total = 0;
	int pass = 0;
	string s1 = "zcat ";
	string s2 = r2 + "_R1.fq.gz";
	string s3 = s1 + s2;
	FILE * red1;
	red1 = popen(s3.c_str(), "r");
	s1 = "gzip - > ";
	s2 = r2 + "_combined.fq.gz";
	s3 = s1 + s2;
	FILE * outfile;
	outfile = popen(s3.c_str(), "w");
	s1 = "zcat ";
	s2 = r2 + "_R2.fq.gz";
	s3 = s1 + s2;
	FILE * red2;
	red2 = popen(s3.c_str(), "r");
	char buffer[2000];
	fqline in_line1;
	fqline in_line2;
	while(fgets(buffer, sizeof(buffer), red1)){
		++total;
		string line1(buffer);
		fgets(buffer, sizeof(buffer), red2);
		string line2(buffer);
		line1 = cxstring::chomp(line1);
		line2 = cxstring::chomp(line2);
		in_line1.read_part_record(red1, line1);
		in_line2.read_part_record(red2, line2);
		read2 read_2;
		read_2.init(in_line2.seq);
		read_2.trim();
		string new_seq = read_2.sbc1 + read_2.sbc2 + read_2.sbc3 + read_2.sbc4;
		string umi = read_2.umi;
		//string umi = in_line.seq.substr(0, 10);
		in_line2.seq = new_seq;
		in_line2.qual = in_line2.qual.substr(0, in_line2.seq.length());
		if(in_line2.seq.length()!=24)continue;
		//in_line2.seq = new_seq;
		string a;
		stringstream as;
		as << in_line2.readname;
		as >> a;
		in_line2.readname = a + ":" + umi + ":" + in_line1.seq + ":" + in_line1.qual;
		//in_line.readname = in_line.readname + ":" + umi;
		in_line2.write_record(outfile);
		++pass;
	}
	pclose(red1);
	pclose(red2);
	pclose(outfile);

	cout << total << " read pairs processed." << endl;
	cout << pass << " read pairs passed docking rate." << endl;
	return;

}
void combine_help(){
	cout << "reachtools combine sample_Prefix" << endl;
}

void pre_process_help(){
	cout << "reachtools pre_process Read2.fq.gz" << endl;
	return;
}


void pairBarcodes(string bam, string target){
	// map <readname, cell_id>
	// map <readname, umi>
	map<string, string> cid;
	map<string, string>::iterator i_cid;

	FILE * infile;
	{
		///if test
		ifstream iftest;
		iftest.open(bam);
		if(!iftest){
			cout << "Cannot open bc_align_file file\n";
			exit(-1);
		}
		iftest.close();
	}
	string s1 = "samtools view ";

	// read bc_align_file to extract bc1-3
	if(bam.substr(bam.length()-3, 3) == "bam"){
		s1 = "samtools view ";

	}
	else if(bam.substr(bam.length()-2, 2) == "gz"){
		s1 = "zcat ";
	}
	string s2 = bam;
	string s3 = s1 + s2;
	infile = popen(s3.c_str(), "r");
	char buffer[2000];
	while(fgets(buffer, sizeof(buffer), infile)){
		string line(buffer);
		line = cxstring::chomp(line);
		samline align;
		align.init(line);
		/*
		string readname = align.readname;
		string cellid = align.chr;
		cid[readname] = cellid;
		*/
		string readname = align.readname.substr(0, align.readname.length()-11);
		string umi = align.readname.substr(align.readname.length()-10, 10);
		cid[readname] = align.chr + ":" + umi;
	}

	// proc target sam/bam file
	{
		///if test
		ifstream iftest;
		iftest.open(target);
		if(!iftest){
			cout << "Cannot open read2 file\n";
			exit(-1);
		}
		iftest.close();
	}
	s1 = "samtools view -h ";
	s2 = target;
	s3 = s1 + s2;
	infile = popen(s3.c_str(), "r");
	FILE * outfile;

	s1 = "samtools view -b - > ";
	s2 = target.substr(0, target.length() - 4) + "_UMI_BC.bam";
	s3 = s1 + s2;
	outfile = popen(s3.c_str(), "w");
	while(fgets(buffer, sizeof(buffer), infile)){
		string line(buffer);
		line = cxstring::chomp(line);
		if(line.substr(0, 1) == "@"){
			fputs((line+"\n").c_str(), outfile);
			continue;
		}

		samline smline;
		smline.init(line);
		i_cid = cid.find(smline.readname);
		if(i_cid == cid.end())continue;
		smline.readname = smline.readname + ":" + i_cid->second;
		smline.write(outfile);
	}
	pclose(outfile);
	pclose(infile);

	return;
}


void pairBarcodes_help(){
	cout << "reachtools pariBarcodes bc_align_file target.sam/bam" << endl;

	return;
}

void pairBC2R1_help(){
	cout << "reachtools pairBarcodes bc_align_file target.sam/bam" << endl;

	return;
}

void convert(string prefix){
	int total = 0;
	int pass = 0;
	string s1 = "samtools view ";
	string s2 = prefix;
	string s3 = s1 + s2;
	FILE * inbam;
	inbam = popen(s3.c_str(), "r");
	s1 = "gzip - > ";
	s2 = prefix.substr(0, prefix.length()-4) + "_cov.fq.gz";
	s3 = s1 + s2;
	FILE * fout;
	fout = popen(s3.c_str(), "w");
	samline align_line;
	fqline fastq_line;
	char buffer[1000];
	while(fgets(buffer, sizeof(buffer), inbam)){
		++total;
		string line(buffer);
		line=cxstring::chomp(line);
		align_line.init(line);
		if(align_line.chr == "*")continue;
		vector<string> tmp = cxstring::split(align_line.readname, ":");
		fastq_line.readname = "@" + tmp[0] + ":" + tmp[1] + ":" + tmp[2] + ":" + tmp[3] + ":" + tmp[4] + ":" + tmp[5] + ":" + tmp[6] + ":" +  align_line.chr + ":" + tmp[7];
		fastq_line.seq = tmp[8];
		fastq_line.qual = tmp[9];
		fastq_line.mark = "+";
		fastq_line.write_record(fout);
		++pass;
	}
	pclose(inbam);
	pclose(fout);
	cout << total << " reads processed." << endl;
	cout << pass << " mapped reads." << endl;
	return;
}
void convert_help(){
	cout << "reachtools convert BC_align.sam" << endl;
	return;

}

