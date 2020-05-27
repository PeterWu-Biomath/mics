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

bool debug = false;

void parse_scaff(char* scaffold_file, std::vector<std::string>& scaff_name, 
	std::vector<std::string>& uncut_scaffs, std::vector<std::vector<std::pair<int, int>>>& scaff_region)
{
	std::ifstream scaff_stream;
	scaff_stream.open(scaffold_file);
	//std::vector < std::string > uncut_scaffs;
	std::string line;
	int count = -1;
	while (!scaff_stream.eof())
	{
		std::getline(scaff_stream, line);
		if (line.find(">") != std::string::npos)
		{
			count++;
			std::string s = line.substr(line.find(">") + 1);

			scaff_name.push_back(s.substr(0, s.find(' ')));
			uncut_scaffs.push_back("");
			//scaffolds.push_back({});
			scaff_region.push_back({});
			//std::cout<<"!"<<std::endl;
		}
		else
		{
			std::remove(line.begin(), line.end(), '\n');
			uncut_scaffs[count].append(line);
		}
	}

	/*if (debug)
	{
		for (auto& scaff : uncut_scaffs)
		{
			std::cout << scaff << std::endl;
		}
	}*/

	for (int i = 0; i <= count; i++)
	{
		auto& scaff = uncut_scaffs[i];
		//scaff.append("N");
		bool n_region = false;
		int first_n = 0;
		for (int j = 0; j < scaff.size(); j++)
		{
			if (!n_region && (scaff[j] == 'n' || scaff[j] == 'N'))
			{
				n_region = true;

				first_n = j;
				//scaff_region[i].push_back(std::make_pair(first_not_n, j));
				//scaffolds[i].push_back(scaff.substr(first_not_n, j - first_not_n));
			}
			else if (n_region && (scaff[j] != 'n' && scaff[j] != 'N'))
			{
				n_region = false;
				scaff_region[i].push_back(std::make_pair(first_n, j));
				//first_not_n = j;
			}
		}
	}

	
}

struct alignment
{
	std::string qname;
	int qlength;
	int qstart;
	int qend;
	bool dir;
	std::string tname;
	int tlength;
	int tstart;
	int tend;
	int nmatch;
	int nmatch_gap;
	int mapq;
	bool prim;
};

void print_align(alignment& a)
{
	std::cout << "alignment:  " << a.qname << " " << a.qlength << " " << a.qstart << " " << a.qend << " "
		<< a.dir << " " << a.tname << " " << a.tstart << " " << a.tend << " " << a.prim << " " << std::endl;
}


void parse_paf(const char* paf_file, std::vector<alignment>& alns)
{
	alignment aln;
	std::ifstream paf_stream;
	paf_stream.open(paf_file);
	std::string line;
	std::string scaff1, scaff2;
	int length1, length2;
	char dir;
	int qstart, qend;
	int tstart, tend;
	int rmatch;
	int align_length;
	int mapq;
	std::string sam;
	std::stringstream linestream;

	while (!paf_stream.eof())
	{
		std::getline(paf_stream, line);
		linestream.str(line);
		linestream >> scaff1 >> length1 >> qstart >> qend >> dir >> scaff2 >> length2 >> tstart >> tend >> rmatch >> align_length >> mapq >> sam;
		aln.qname = scaff1;
		aln.qlength = length1;
		aln.qstart = qstart;
		aln.qend = qend;
		aln.dir = (dir == '+');
		aln.tname = scaff2;
		aln.tlength = length2;
		aln.tstart = tstart;
		aln.tend = tend;
		aln.nmatch = rmatch;
		if (rmatch < 200) continue;
		aln.nmatch_gap = align_length;
		aln.mapq = mapq;
		aln.prim = (line.find("tp:A:P") != std::string::npos);
		alns.push_back(aln);
	}
	//for (auto& aln : alns) { print_align(aln); }
}
void printseq(std::string s, int step)
{
	//std::cout << s << std::endl;
	for (int i = 0; i < s.size(); i += step)
	{
		if (i + step >= s.size()) {
			std::cout << s.substr(i, s.size()-i) << std::endl;
		}
		else {
			std::cout << s.substr(i, step) << std::endl;
		}
	}
}

int main(int argc, char** argv)
{
	char* scaffold_file = argv[1];

	//std::vector<std::vector<std::string>> scaffolds;
	std::vector<std::string> scaff_name,scaffs;
	std::vector<std::vector<std::pair<int, int>>> scaff_region;

	parse_scaff(scaffold_file, scaff_name,scaffs, scaff_region);

	/*for (int i = 0; i < scaff_name.size(); i++)
	{
		for (int j = 0; j < scaff_region[i].size(); j++)
		{
			std::cout << scaff_name[i] << "\t" << scaff_region[i][j].first 
				<< "\t" << scaff_region[i][j].second << std::endl;
		}
	}*/
	std::vector<alignment> alns;
	parse_paf(argv[2], alns);

	std::map<std::string, int> scaff_name_inv;
	for (int i = 0; i < scaff_name.size(); i++) { scaff_name_inv[scaff_name[i]] = i; }

	int cnt = 0;
	for (int i = 1; i < alns.size(); i++)
	{
		auto& aln1 = alns[i - 1];
		auto& aln2 = alns[i];
		if (aln1.qname != aln2.qname) continue;
		//if (aln1.qname != aln2.qname) continue;
		//if (aln1.dir != aln2.dir) continue;
		if (aln1.qend > aln2.qstart - 300) continue;
		
		int id = scaff_name_inv[aln2.qname];
		if (id < 0) {
			std::cerr <<"?"<< aln2.qname << std::endl; continue;
		}
		auto& nregion = scaff_region[id];
		bool intersect=false;
		for (int j=0;j< nregion.size() &&!intersect;j++)
		{
			auto& p = nregion[j];
			if (aln1.qend < p.second && p.first < aln2.qstart) { intersect = true; break; }
		}
		if (!intersect)
		{
			/*if (aln1.dir)
			{
				//if (aln2.qstart - aln1.qend < 200)
				{
					print_align(aln1);
					print_align(aln2);
					std::cout << std::endl;
				}
			}
			else //if(aln2.qstart - aln1.qend<200)
			{
				print_align(aln1);
				print_align(aln2);
				std::cout << std::endl;
			}*/
			cnt++;
			std::cout << ">" << cnt << " " << aln1.qname << " " << aln1.qend << "->" << aln2.qstart
				<<" "<<aln1.tname<<"@"<< (aln1.dir? aln1.tend:aln1.tstart)<<"->"
				<< aln2.tname << "@" << (aln2.dir ? aln2.tstart : aln2.tend) << std::endl;
			printseq(scaffs[id].substr(aln1.qend, aln2.qstart- aln1.qend), 80);
		}
	}

	return 0;
}
