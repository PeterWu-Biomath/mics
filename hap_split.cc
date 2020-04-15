#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <utility>
#include <climits>
#include <omp.h>
#include <tuple>

#include "htslib/sam.h"

#include "ozstream.hpp"

using namespace std;
using namespace zstream;

std::vector<int> base2num(256, -1);
std::vector<char> num2base = { 'A','C','G','T' };

std::vector<int> base2numx(16, -1);

std::string Date()
{
	time_t now = time(0);

	// convert now to string form
	std::string dt = (std::string)(ctime(&now));
	dt.resize(dt.size() - 1);
	return dt + ": ";
}

struct snp
{
	int chr;
	int pos;
	std::vector<int> ref_base;
	std::vector<int> alt_base;
	snp(int chr_, int pos_, std::vector<int> ref_base_, std::vector<int> alt_base_)
	{
		chr = chr_;
		pos = pos_;
		ref_base = ref_base_;
		alt_base = alt_base_;
	}
	const bool operator < (const snp& b) const {
		if (chr != b.chr)
		{
			return chr < b.chr;
		}
		else if (pos != b.pos)
		{
			return pos < b.pos;
		}
		else if (ref_base != b.ref_base)
		{
			return ref_base < b.ref_base;
		}
		else { return alt_base < b.alt_base; }
	}

	const bool operator == (const snp& b) const {
		return ((chr == b.chr) && (pos == b.pos) && (ref_base == b.ref_base) && (alt_base == b.alt_base));
	}
	const bool operator > (const snp& b)
	{
		return (b < *this);
	}
};

struct snp_diff
{
	int chr;
	int pos;
	std::vector<int> m_base;
	std::vector<int> p_base;
	snp_diff(int chr_, int pos_, std::vector<int> m_base_, std::vector<int> p_base_)
	{
		chr = chr_;
		pos = pos_;
		m_base = m_base_;
		p_base = p_base_;
	}
	const bool operator < (const snp_diff& b) const {
		if (chr != b.chr)
		{
			return chr < b.chr;
		}
		else if (pos != b.pos)
		{
			return pos < b.pos;
		}
		else if (m_base != b.m_base)
		{
			return m_base < b.m_base;
		}
		else { return p_base < b.p_base; }
	}

	const bool operator == (const snp_diff& b) const {
		return ((chr == b.chr) && (pos == b.pos) && (m_base == b.m_base) && (p_base == b.p_base));
	}
	const bool operator > (const snp_diff& b)
	{
		return (b < *this);
	}
};


struct {
	bool operator()(snp a, snp b) const
	{
		if (a.chr != b.chr)
		{
			return a.chr < b.chr;
		}
		else if (a.pos != b.pos)
		{
			return a.pos < b.pos;
		}
		else if (a.ref_base != b.ref_base)
		{
			return a.ref_base < b.ref_base;
		}
		else { return a.alt_base < b.alt_base; }
	}
} snpLess;

std::vector<int> parsebases(std::string seq)
{
	std::vector<int> v;
	for (int i = seq.size() - 1; i >= 0; i--)
	{
		v.push_back(base2num[seq[i]]);

		/*if (seq[i] == 'A' || seq[i] == 'a') { v.push_back(0); }
		else if (seq[i] == 'T' || seq[i] == 't') { v.push_back(3); }
		else if (seq[i] == 'C' || seq[i] == 'c') { v.push_back(1); }
		else if (seq[i] == 'G' || seq[i] == 'g') { v.push_back(2); }*/
	}
	return v;
}

void Read_snp(const char* file, std::vector <snp>& snps)
{
	std::string line;
	std::ifstream myfile;
	myfile.open(file);

	//(std::string(file));
	//if (myfile.is_open())
	//{
	while (!myfile.eof())
	{
		std::getline(myfile, line);
		if (line[0] > '9' || line[0] < '1') continue;
		int chr, pos;
		std::string ref, alt, comment;
		std::stringstream s(line); // Used for breaking words 
		s >> chr >> pos >> comment >> ref >> alt;
		if (ref.size() != 1) continue;
		if (alt.size() != 1) continue;
		snps.push_back(snp(chr - 1, pos - 1, parsebases(ref), parsebases(alt)));

		/*std::cout << chr << "\t"
			<< pos << "\t"
			<< ref << "\t"
			<< alt << std::endl;*/
	}
	myfile.close();
	std::sort(snps.begin(), snps.end());
	//}
}

bool Substr(std::string a, std::string b)
{
	return (a.find(b) != std::string::npos);
}

int Het_comment(std::string a)
{
	if (Substr(a, "0|1")) return 1;
	if (Substr(a, "1|0")) return 2;
	if (Substr(a, "2|1")) return 4;
	if (Substr(a, "1|2")) return 3;
	return -1;
}

void Read_snp_diff(const char* file,
	const std::vector< std::string>& scaffolds,
	std::vector <snp_diff>& snp_diffs)
{
	snp_diffs.reserve(1000000);
	std::string line;
	std::ifstream myfile;
	myfile.open(file);

	while (!myfile.eof())
	{
		std::getline(myfile, line);
		if (line[0] > '9' || line[0] < '1') continue;
		int chr, pos;
		std::string ref, alt, id, qual, filter, info, format, comment;
		std::stringstream s(line); // Used for breaking words 
		s >> chr >> pos >> id >> ref >> alt >> qual >> filter >> info >> format >> comment;
		if (ref.size() != 1) continue;
		if (alt.size() != 1 && alt.size() != 3) continue;
		int het = Het_comment(line);
		if (het==-1) continue;
		if (alt.size() == 1)
		{
			if (het == 1)
			{
				snp_diffs.push_back(snp_diff(chr - 1, pos - 1, 
					{ base2num[ref[0]] }, { base2num[alt[0]] }));
			}
			else
			{
				snp_diffs.push_back(snp_diff(chr - 1, pos - 1,
					{ base2num[alt[0]] }, { base2num[ref[0]] }));
			}
		}
		if (alt.size() == 3)
		{
			std::vector<int> a = { base2num[alt[0]] };
			std::vector<int> b = {base2num[alt[3]]};
			if (het == 3)
			{
				snp_diffs.push_back(snp_diff(chr - 1, pos - 1,
					a,b));
			}
			else
			{
				snp_diffs.push_back(snp_diff(chr - 1, pos - 1,
					b,a));
			}
		}
		//snps.push_back(snp(chr - 1, pos - 1, parsebases(ref), parsebases(alt)));
	}

	myfile.close();
	std::sort(snp_diffs.begin(), snp_diffs.end());
}

void parse_scaff(const char* scaffold_file,
	std::vector< std::string>& scaffolds)
{
	std::vector<std::string> scaff_name;
	std::ifstream scaff_stream;
	scaff_stream.open(scaffold_file);
	if (!scaff_stream)
	{
		std::cerr << "Open Scaffold file failed !" << std::endl;
	}
	std::string line;
	//std::string tmp_scaff_name;
	int count = -1;
	while (!scaff_stream.eof())
	{
		std::getline(scaff_stream, line);
		if (line.find(">") != std::string::npos)
		{
			std::remove(line.begin(), line.end(), '\n');
			scaff_name.push_back(line.substr(line.find(">") + 1));
			scaffolds.push_back("");
			scaffolds.back().reserve(100000000);
			count++;
		}
		else
		{
			std::remove(line.begin(), line.end(), '\n');
			scaffolds[count].append(line);
		}
	}
}

void bam_get_read_seq(bam1_t* b, std::string& str)
{
	int i = 0;
	uint8_t* s = bam_get_seq(b);
	for (; i < b->core.l_qseq; ++i)
	{
		str += "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
	}

}
std::string ToString(bam1_t* bam_record)
{
	std::string a;
	a += "@";
	a += bam_get_qname(bam_record);
	a += "\n";
	uint8_t* s = bam_get_seq(bam_record);
	for (int i = 0; i < bam_record->core.l_qseq; ++i)
	{
		a += "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
	}
	a += "\n+\n";
	for (int i = 0; i < bam_record->core.l_qseq; ++i)
	{
		a += (char)(bam_get_qual(bam_record)[i] + 33);
	}
	a += "\n";
	return a;
}

void StatErrorProfile(char* bamfile,
	const std::vector< std::string>& ref,
	std::vector<std::vector<double>>& mismatch,
	std::vector<double>& ins,
	std::vector<double>& del)
{
	std::vector<std::vector<int64_t>> mismatch_cnt(4, { 0,0,0,0 });
	std::vector<int64_t> ins_cnt(4, 0);
	std::vector<int64_t> del_cnt(4, 0);
	int64_t total;
	samFile* bam_fp = nullptr;
	bam1_t* bam_record = bam_init1();

	mismatch = { { 0,0,0,0 } ,{ 0,0,0,0 } ,{ 0,0,0,0 } ,{ 0,0,0,0 } };
	ins = { 0,0,0,0 };
	del = { 0,0,0,0 };

	//std::vector<bam1_t> alns;

	// open bam file.0
	if ((bam_fp = sam_open(bamfile, "r")) == nullptr) {
		std::cerr << "Can't open bam file: " << bamfile << std::endl;
	}

	// get header info.
	bam_hdr_t* bam_header = bam_hdr_init();
	if ((bam_header = sam_hdr_read(bam_fp)) == nullptr) {
		std::cerr << "Can't read bam header." << std::endl;
	}

	int64_t read_cnt = 0, read_filtered = 0, cnt = 0;
	while (sam_read1(bam_fp, bam_header, bam_record) >= 0 && cnt < 100000)
	{
		cnt++;
		if (cnt % 10000 == 0) { std::cout << Date() << "Processed " << cnt << " reads." << std::endl; }
		bam1_t& aln = *bam_record;
		//std::cout << ToString(bam_record);
		if ((aln.core.flag & BAM_FUNMAP != 0) || (aln.core.flag & BAM_FSECONDARY != 0)) continue;
		int rstart = aln.core.pos;
		int qstart = 0;
		bool inv = !bam_is_rev(&aln);
		auto& ref_seq = ref[aln.core.tid];
		for (int i = 0; i < aln.core.n_cigar; i++)
		{
			auto c = bam_get_cigar(&aln)[i];
			uint32_t l = c / 16;

			if (c % 16 == 0 || c % 16 == 8 || c % 16 == 7)
			{
				for (int loop = 0; loop < l; loop++)
				{
					auto& ref_base = ref_seq[rstart + loop];
					int q_base = bam_seqi(bam_get_seq(&aln), qstart + loop);

					int ind1 = base2num[ref_base], ind2 = base2numx[q_base];
					if (ind1 >= 0 && ind2 >= 0)
					{
						if (inv) { ind1 = 3 - ind1; ind2 = 3 - ind2; }
						mismatch_cnt[ind1][ind2]++;
						total++;
					}
					/*if (ref_base == 'A' || ref_base == 'a') { ind1 = 0; }
					else if (ref_base == 'C' || ref_base == 'c') { ind1 = 1; }
					else if (ref_base == 'G' || ref_base == 'g') { ind1 = 2; }
					else if (ref_base == 'T' || ref_base == 't') { ind1 = 3; }
					if (q_base == 1) {}*/

				}

				rstart += l;
				qstart += l;
			}

			else if (c % 16 == 4) { qstart += l; }
			else if (c % 16 == 1)
			{
				for (int loop = 0; loop < l; loop++)
				{
					int q_base = bam_seqi(bam_get_seq(&aln), qstart + loop);

					int  ind2 = base2numx[q_base];
					if (ind2 >= 0)
					{

						if (inv) { ind2 = 3 - ind2; }
						ins_cnt[ind2]++;
						total++;
					}

				}
				qstart += l;
			}

			else if (c % 16 == 2)
			{
				for (int loop = 0; loop < l; loop++)
				{
					auto& ref_base = ref_seq[rstart + loop];
					int ind1 = base2num[ref_base];
					if (ind1 >= 0)
					{
						if (inv) { ind1 = 3 - ind1; }
						del_cnt[ind1]++;
						total++;
					}
				}
				rstart += l;
			}
		}

		/*if (bam_is_rev(&aln)) { std::cout << "RC\t"; }
		else { std::cout << "+\t"; }
		std::string base;
		bam_get_read_seq(&aln, base);
		std::cout << base << std::endl;*/
	}
	sam_close(bam_fp);

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mismatch[i][j] = double(mismatch_cnt[i][j]) /
				(mismatch_cnt[i][0] + mismatch_cnt[i][1]
					+ mismatch_cnt[i][2] + mismatch_cnt[i][3]);

			std::cout << mismatch[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	for (int i = 0; i < 4; i++)
	{

		ins[i] = double(ins_cnt[i]) / total;
		del[i] = double(del_cnt[i]) /
			(mismatch_cnt[i][0] + mismatch_cnt[i][1]
				+ mismatch_cnt[i][2] + mismatch_cnt[i][3]);

		std::cout << ins[i] << "\t" << del[i] << std::endl;
	}

}
void init()
{
	base2num['A'] = 0;
	base2num['a'] = 0;
	base2num['C'] = 1;
	base2num['c'] = 1;
	base2num['G'] = 2;
	base2num['g'] = 2;
	base2num['T'] = 3;
	base2num['t'] = 3;


	base2numx[1] = 0;
	base2numx[2] = 1;
	base2numx[4] = 2;
	base2numx[8] = 3;
}

void Splitbam(char* bamfile, const std::string prefix,
	const std::vector< std::string>& ref,
	const std::vector<snp_diff>& snp_diff_i,
	const std::vector<std::vector<double>>& mismatch,
	const std::vector<double>& ins,
	const std::vector<double>& del)
{
	int chr_num = snp_diff_i.back().chr + 1;
	std::vector< std::vector<uint64_t>> ind(chr_num);

	{
		for (int chr = 0; chr < chr_num; chr++)
		{
			//std::cout << Date() << ref[chr].size() << std::endl;
			ind[chr].resize(ref[chr].size(), 0);

			//std::cout << Date() <<"Finish "<<ref[chr].size() << std::endl;
			//ind2[chr].resize(ref[chr], 0);
		}
		std::cout << Date() << "Reserving space complete" << std::endl;

		size_t start = 0;
		for (size_t chr = 0; chr < chr_num; chr++)
		{
			auto& ind_i = ind[chr];
			for (size_t loop = 0; loop < snp_diff_i[start].pos; loop++)
			{
				ind_i[loop] = start;
			}

			for (size_t loop = start + 1; loop < snp_diff_i.size(); loop++)
			{
				if (snp_diff_i[loop].chr != chr) {
					start = loop; break;
				}
				for (size_t loop1 = snp_diff_i[loop - 1].pos; loop1 < snp_diff_i[loop].pos; loop1++)
				{
					ind_i[loop1] = loop;
				}
			}
			for (size_t loop = snp_diff_i[start - 1].pos; loop < ind_i.size(); loop++)
			{
				ind_i[loop] = start;
			}
			//std::cout << Date() << "Build index for chr"<<chr<<" complete." << std::endl;

		}
	}
	std::cout << Date() << "Building index complete" << std::endl;
	std::cout << Date() << "Start phaing." << std::endl;

	std::ofstream out1x, out2x;
	out1x.open(prefix + ".m.fastq.gz");
	out2x.open(prefix + ".p.fastq.gz");
	//ofstream os("test.txt.gz");
	ogzstream out1(out1x);
	ogzstream out2(out2x);
	//reading bam
	samFile* bam_fp = nullptr;
	bam1_t* bam_record = bam_init1();

	//std::vector<bam1_t> alns(8192);
	std::vector<std::tuple<int, int, int>> score;
	score.reserve(1000000);

	// open bam file.0
	if ((bam_fp = sam_open(bamfile, "r")) == nullptr) {
		std::cerr << "Can't open bam file: " << bamfile << std::endl;
	}

	// get header info.
	bam_hdr_t* bam_header = bam_hdr_init();
	if ((bam_header = sam_hdr_read(bam_fp)) == nullptr) {
		std::cerr << "Can't read bam header." << std::endl;
	}

	int64_t cnt = 0;
	while (sam_read1(bam_fp, bam_header, bam_record) >= 0 )
	{
		/*alns[cnt % 8192] = *bam_record;
		if (cnt % 8192 == 8191)
		{
//#pragma omp parallel for num_threads(32)
			for (int i = 0; i < 8192; i++)
			{*/
			//std::cout << "start" << std::endl;
			//score[i] = std::make_tuple(0, 0, 0);
		auto& aln = *bam_record;
		if ((aln.core.flag & 256) != 0 || (aln.core.flag & 2048) != 0)
		{
			//score[i] = std::make_tuple(0, 0, -1);
			continue;
		}
		if ((aln.core.flag & 4) != 0)
		{
			out1 << ToString(bam_record);
			out2 << ToString(bam_record);
			continue;
		}
		if (aln.core.tid >= chr_num)
		{
			if (aln.core.tid == chr_num)
			{
				out1 << ToString(bam_record);
			}
			else if (aln.core.tid == chr_num)
			{
				out2 << ToString(bam_record);
			}
			else
			{
				out1 << ToString(bam_record);
				out2 << ToString(bam_record);
			}
			continue;
		}
		//size_t snp_start = ind[aln.core.tid][aln.core.pos < 1 ? 0 : (aln.core.pos - 1)];

		int chr = aln.core.tid;
		int rpos = aln.core.pos, qpos = 0;

		int total = 0, m_cnt = 0, p_cnt = 0;
		for (int loop = 0; loop < aln.core.n_cigar; loop++)
		{
			auto c = bam_get_cigar(&aln)[loop];

			uint32_t l = c / 16;
			//c = c % 16;
			/*std::cout << c << "\t"
				<< l << "\t"
				<< c % 16 << std::endl;*/
			if (c % 16 == 0 || c % 16 == 8 || c % 16 == 7)
			{
				size_t snp_start = ind[aln.core.tid][rpos < 1 ? 0 : (rpos - 1)];

				size_t snp_stop = ind[aln.core.tid][rpos + l - 1];

				for (int j = snp_start; j < snp_stop; j++)
				{
					total++;
					int m = snp_diff_i[j].m_base[0];
					int p = snp_diff_i[j].p_base[0];
					if (snp_diff_i[j].pos - rpos + qpos > aln.core.l_qseq) continue;
					int base = base2numx[bam_seqi(bam_get_seq(&aln), snp_diff_i[j].pos - rpos + qpos)];
					if (base == m) { m_cnt++; }
					if (base == p) { p_cnt++; }
				}

				rpos += l;
				qpos += l;
			}

			else if (c % 16 == 4) { qpos += l; }
			else if (c % 16 == 1) { qpos += l; }
			else if (c % 16 == 2) { rpos += l; }

		}
		std::cout << "Score:\t" << m_cnt << "\t" << p_cnt << "\t" << total << std::endl;
		if (m_cnt > p_cnt * 3 + 1) {
			out1 << ToString(bam_record);
		}
		else if (p_cnt > m_cnt * 3 + 1) {
			out2 << ToString(bam_record);
		}
		else
		{

			out1 << ToString(bam_record);
			out2 << ToString(bam_record);
		}
		if (cnt % 1000 == 0)
		{
			std::cout << Date() << "Splited " << cnt << " reads." << std::endl;
		}
		cnt++;
	}
	sam_close(bam_fp);

	out1.close();

	out2.close();
}

void Splitbam(char* bamfile, const std::string prefix,
	const std::vector< std::string>& ref,
	const std::vector<snp>& snp_m,
	const std::vector<snp>& snp_p,
	const std::vector<std::vector<double>>& mismatch,
	const std::vector<double>& ins,
	const std::vector<double>& del)
{
	int chr_num = snp_m.back().chr + 1;
	std::vector< std::vector<uint64_t>> ind(chr_num);
	std::vector<snp_diff> snp_diff_i;
	snp_diff_i.reserve(snp_m.size() + snp_p.size());
	//build diff
	{
		size_t i1 = 0, i2 = 0;
		while (i1 < snp_m.size() && i2 < snp_p.size())
		{
			if ((snp_m[i1].chr < snp_p[i2].chr)
				|| ((snp_m[i1].chr == snp_p[i2].chr) && (snp_m[i1].pos < snp_p[i2].pos)))
			{
				snp_diff_i.push_back(snp_diff(snp_m[i1].chr,
					snp_m[i1].pos, snp_m[i1].alt_base, snp_m[i1].ref_base));
				i1++;
			}
			else if ((snp_p[i2].chr < snp_m[i1].chr)
				|| ((snp_m[i1].chr == snp_p[i2].chr) && (snp_m[i1].pos > snp_p[i2].pos)))
			{
				snp_diff_i.push_back(snp_diff(snp_p[i2].chr,
					snp_p[i2].pos, snp_p[i2].ref_base, snp_p[i2].alt_base));
				i2++;
			}
			else
			{
				if (snp_m[i1].alt_base != snp_p[i2].alt_base)
				{
					snp_diff_i.push_back(snp_diff(snp_p[i2].chr,
						snp_p[i2].pos, snp_m[i1].alt_base, snp_p[i2].alt_base));
				}
				i1++;
				i2++;
			}
		}

		for (; i1 < snp_m.size(); i1++)
		{
			snp_diff_i.push_back(snp_diff(snp_m[i1].chr,
				snp_m[i1].pos, snp_m[i1].alt_base, snp_m[i1].ref_base));
		}
		for (; i2 < snp_p.size(); i2++)
		{
			snp_diff_i.push_back(snp_diff(snp_p[i2].chr,
				snp_p[i2].pos, snp_p[i2].ref_base, snp_p[i2].alt_base));
		}
	}
	std::cout << Date() << "Building diff complete" << std::endl;

	/*for (auto& diff : snp_diff_i)
	{
		std::cout << diff.chr << "\t"
			<< diff.pos << "\t"
			<< num2base[diff.m_base[0]] << "\t"
			<< num2base[diff.p_base[0]] << std::endl;
	}*/

	//build index 
	{
		for (int chr = 0; chr < chr_num; chr++)
		{
			//std::cout << Date() << ref[chr].size() << std::endl;
			ind[chr].resize(ref[chr].size(), 0);

			//std::cout << Date() <<"Finish "<<ref[chr].size() << std::endl;
			//ind2[chr].resize(ref[chr], 0);
		}
		std::cout << Date() << "Reserving space complete" << std::endl;

		size_t start = 0;
		for (size_t chr = 0; chr < chr_num; chr++)
		{
			auto& ind_i = ind[chr];
			for (size_t loop = 0; loop < snp_diff_i[start].pos; loop++)
			{
				ind_i[loop] = start;
			}

			for (size_t loop = start + 1; loop < snp_diff_i.size(); loop++)
			{
				if (snp_diff_i[loop].chr != chr) {
					start = loop; break;
				}
				for (size_t loop1 = snp_diff_i[loop - 1].pos; loop1 < snp_diff_i[loop].pos; loop1++)
				{
					ind_i[loop1] = loop;
				}
			}
			for (size_t loop = snp_diff_i[start - 1].pos; loop < ind_i.size(); loop++)
			{
				ind_i[loop] = start;
			}
			//std::cout << Date() << "Build index for chr"<<chr<<" complete." << std::endl;

		}
	}
	std::cout << Date() << "Building index complete" << std::endl;
	std::cout << Date() << "Start phaing." << std::endl;

	std::ofstream out1, out2;
	out1.open(prefix + ".m.fastq");
	out2.open(prefix + ".p.fastq");

	//reading bam
	samFile* bam_fp = nullptr;
	bam1_t* bam_record = bam_init1();

	//std::vector<bam1_t> alns(8192);
	std::vector<std::tuple<int, int, int>> score;
	score.reserve(1000000);

	// open bam file.0
	if ((bam_fp = sam_open(bamfile, "r")) == nullptr) {
		std::cerr << "Can't open bam file: " << bamfile << std::endl;
	}

	// get header info.
	bam_hdr_t* bam_header = bam_hdr_init();
	if ((bam_header = sam_hdr_read(bam_fp)) == nullptr) {
		std::cerr << "Can't read bam header." << std::endl;
	}

	int64_t cnt = 0;
	while (sam_read1(bam_fp, bam_header, bam_record) >= 0)
	{
		/*alns[cnt % 8192] = *bam_record;
		if (cnt % 8192 == 8191)
		{
//#pragma omp parallel for num_threads(32)
			for (int i = 0; i < 8192; i++)
			{*/
			//std::cout << "start" << std::endl;
			//score[i] = std::make_tuple(0, 0, 0);
		auto& aln = *bam_record;
		if ((aln.core.flag & 256) != 0 || (aln.core.flag & 2048) != 0)
		{
			//score[i] = std::make_tuple(0, 0, -1);
			continue;
		}
		if (aln.core.flag & 4 != 0)
		{
			out1 << ToString(bam_record);
			out2 << ToString(bam_record);
			continue;
		}
		if (aln.core.tid >= chr_num)
		{
			if (aln.core.tid == chr_num)
			{
				out1 << ToString(bam_record);
			}
			else if (aln.core.tid == chr_num)
			{
				out2 << ToString(bam_record);
			}
			else
			{
				out1 << ToString(bam_record);
				out2 << ToString(bam_record);
			}
			continue;
		}
		//size_t snp_start = ind[aln.core.tid][aln.core.pos < 1 ? 0 : (aln.core.pos - 1)];

		int chr = aln.core.tid;
		int rpos = aln.core.pos, qpos = 0;

		int total = 0, m_cnt = 0, p_cnt = 0;
		for (int loop = 0; loop < aln.core.n_cigar; loop++)
		{
			auto c = bam_get_cigar(&aln)[loop];

			uint32_t l = c / 16;
			//c = c % 16;
			/*std::cout << c << "\t"
				<< l << "\t"
				<< c % 16 << std::endl;*/
			if (c % 16 == 0 || c % 16 == 8 || c % 16 == 7)
			{
				size_t snp_start = ind[aln.core.tid][rpos < 1 ? 0 : (rpos - 1)];

				size_t snp_stop = ind[aln.core.tid][rpos + l - 1];

				for (int j = snp_start; j < snp_stop; j++)
				{
					total++;
					int m = snp_diff_i[j].m_base[0];
					int p = snp_diff_i[j].p_base[0];
					if (snp_diff_i[j].pos - rpos + qpos > aln.core.l_qseq) continue;
					int base = base2numx[bam_seqi(bam_get_seq(&aln), snp_diff_i[j].pos - rpos + qpos)];
					if (base == m) { m_cnt++; }
					if (base == p) { p_cnt++; }
				}

				rpos += l;
				qpos += l;
			}

			else if (c % 16 == 4) { qpos += l; }
			else if (c % 16 == 1) { qpos += l; }
			else if (c % 16 == 2) { rpos += l; }

		}
		std::cout << "Score:\t" << m_cnt << "\t" << p_cnt << "\t" << total << std::endl;
		if (m_cnt > p_cnt * 3 + 1) {
			out1 << ToString(bam_record);
		}
		else if (p_cnt > m_cnt * 3 + 1) {
			out2 << ToString(bam_record);
		}
		else
		{

			out1 << ToString(bam_record);
			out2 << ToString(bam_record);
		}
		if (cnt % 1000 == 0)
		{
			std::cout << Date() << "Splited " << cnt << " reads." << std::endl;
		}
		cnt++;
	}
	sam_close(bam_fp);

	out1.close();

	out2.close();
}

int main(int argc, char** argv)
{

	init();
	//std::vector<snp> snp_m, snp_p;
	//snp_m.reserve(1000000);
	//snp_p.reserve(1000000);
	std::cout << Date() << "Reading Reference" << std::endl;
	std::vector< std::string> scaffolds;
	parse_scaff(argv[2], scaffolds);
	std::cout << Date() << "Reading vcf" << std::endl;
	std::vector<snp_diff> snp_diff_i;
	Read_snp_diff(argv[1], scaffolds, snp_diff_i);
	std::cout << Date() << "Captured "<< snp_diff_i .size()<< " hetero snp." << std::endl;

	std::cout << Date() << "Stating Error Profile." << std::endl;
	std::vector<std::vector<double>> mismatch;
	std::vector<double> ins;
	std::vector<double> del;
	StatErrorProfile(argv[3], scaffolds, mismatch, ins, del);
	Splitbam(argv[3], argv[4], scaffolds, snp_diff_i, mismatch, ins, del);
	
	return 0;
}
