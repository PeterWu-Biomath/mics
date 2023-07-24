#include "htslib/sam.h"
#include "htslib/hts.h"

struct base_diff
{
    int chr;
    hts_pos_t pos;
    std::string ref;
    std::string alt;
};

//a code snippet for parsing cigar and MD 

std::string fetch_base_by_ref_pos(const bam1_t &bam_record, hts_pos_t start, hts_pos_t len)
{
    std::string str(len, '-');
    int cigar_read_pos = 0, cigar_ref_pos = 0;
    uint32_t *cigar = bam_get_cigar(&bam_record);
    for (size_t cigar_id = 0; cigar_id < bam_record.core.n_cigar; cigar_id++)
    {
        if (bam_cigar_type(bam_cigar_op(cigar[cigar_id])) & 1)
        {
            cigar_read_pos += bam_cigar_oplen(cigar[cigar_id]);
        }
        if (bam_cigar_type(bam_cigar_op(cigar[cigar_id])) & 2)
        {
            cigar_ref_pos += bam_cigar_oplen(cigar[cigar_id]);
        }
        // std::cout<<"debug "<<cigar_ref_pos<<' '<<start<<' '<<len<<std::endl;
        if (cigar_ref_pos >= start & bam_cigar_op(cigar[cigar_id]) == 0)
        {
            int start_read_pos = start + cigar_read_pos - cigar_ref_pos;
            for (size_t i = 0; i < len; i++)
            {
                str[i] = seq_nt16_str[bam_seqi(bam_get_seq(&bam_record), start_read_pos + i)];
            }
            break;
        }
    }
    return str;
}

std::vector<base_diff> get_single_diff2(const bam1_t &bam_record)
{
    char *MD_str = (char *)bam_aux_get(&bam_record, "MD") + 1;
    std::vector<base_diff> diff_list;
    int chr_id = bam_record.core.tid;
    hts_pos_t pos = bam_record.core.pos;
    int md_len = 0;
    std::vector<std::pair<int, char *>> MD_list;
    while (*MD_str != '\0')
    {
        if (*MD_str >= '0' & *MD_str <= '9')
        {
            md_len *= 10;
            md_len += (*MD_str - '0');
            MD_str++;
        }
        else
        {
            MD_list.push_back(std::make_pair(md_len, MD_str));
            md_len = 0;
            while (*MD_str<'0' | *MD_str> '9')
            {
                MD_str++;
            }
        }
    }
    MD_list.push_back(std::make_pair(md_len, MD_str));
    int cigar_read_pos = 0, cigar_ref_pos = 0;
    uint32_t *cigar = bam_get_cigar(&bam_record);
    for (size_t cigar_id = 0; cigar_id < bam_record.core.n_cigar; cigar_id++) // consider only M,I,D,S case, other to be implemented
    {
        if (bam_cigar_op(cigar[cigar_id]) == 4)
        {
            cigar_read_pos += bam_cigar_oplen(cigar[cigar_id]);
            // MD_offset += bam_cigar_oplen(cigar[cigar_id]);
        }
        else if (bam_cigar_op(cigar[cigar_id]) == 1) // insertion, get base from reads
        {
            if (cigar_read_pos >= 1) // cigar like 2I94M could happen, no good solution
            {
                std::string alt(bam_cigar_oplen(cigar[cigar_id]) + 1, ' ');
                for (size_t id = 0; id <= bam_cigar_oplen(cigar[cigar_id]); id++)
                {
                    alt[id] = seq_nt16_str[bam_seqi(bam_get_seq(&bam_record), id + cigar_read_pos - 1)];
                }
                std::string ref = std::string(1, seq_nt16_str[bam_seqi(bam_get_seq(&bam_record), cigar_read_pos - 1)]);
                diff_list.push_back({chr_id, pos + cigar_ref_pos, ref, alt});
            }
            cigar_read_pos += bam_cigar_oplen(cigar[cigar_id]);
        }
        else if (bam_cigar_op(cigar[cigar_id]) == 0) // match
        {
            cigar_read_pos += bam_cigar_oplen(cigar[cigar_id]);
            cigar_ref_pos += bam_cigar_oplen(cigar[cigar_id]);
        }
        else if (bam_cigar_op(cigar[cigar_id]) == 2) // del
        {
            if (cigar_read_pos >= 1)
            {
                int len = 0;
                std::string ref, alt;
                ref.append(1, seq_nt16_str[bam_seqi(bam_get_seq(&bam_record), cigar_read_pos - 1)]);
                alt = ref;
                for (size_t id = 0; id < MD_list.size(); id++)
                {
                    len += MD_list[id].first;
                    char *ptr = MD_list[id].second;
                    if (len >= cigar_ref_pos)
                    {
                        // if(len!=cigar_ref_pos){std::cout<<"error"<<std::endl;}
                        ptr++;
                        while (*ptr > '9' & *ptr != '\0')
                        {
                            ref.append(1, *ptr);
                            ptr++;
                        }
                        break;
                    }
                    if (*ptr == '^')
                    {
                        ptr++;
                    }
                    while (*ptr > '9' & *ptr != '\0')
                    {
                        len++;
                        ptr++;
                    }
                }
                diff_list.push_back({chr_id, pos + cigar_ref_pos, ref, alt});
            }
            cigar_ref_pos += bam_cigar_oplen(cigar[cigar_id]);
        }
    }
    // second round, get mismatch from MD_tag.
    hts_pos_t MD_ref_pos = 0;
    for (size_t MD_id = 0; MD_id < MD_list.size(); MD_id++)
    {
        MD_ref_pos += MD_list[MD_id].first;
        hts_pos_t MD_len = 0;
        char *ptr = MD_list[MD_id].second;
        std::string ref, alt;
        if (*ptr == '\0')
        {
            continue;
        }
        else if (*ptr == '^')
        {
            ptr++;
            while (*ptr > '9' & *ptr != '\0')
            {
                MD_ref_pos++;
                ptr++;
            }
            continue;
        }
        while (*ptr > '9' & *ptr != '\0')
        {
            ref.append(1, *ptr);
            MD_len++;
            ptr++;
        }
        // fetch base with ref pos at MD_ref_pos ~ MD_ref_pos+MD_len
        alt = fetch_base_by_ref_pos(bam_record, MD_ref_pos, MD_len);
        diff_list.push_back({chr_id, pos + MD_ref_pos, ref, alt});
        MD_ref_pos += MD_len;
    }
    return diff_list;
}

