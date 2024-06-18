from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from tqdm import tqdm
import os
 
def blast_fasta(fasta_path, db_path, output_path, max_hits=5):
    # 读取输入的FASTA序列
    fasta_sequences = list(SeqIO.parse(fasta_path, "fasta"))
    num_sequences = len(fasta_sequences)
 
    with open(output_path, 'w') as output_file:
        for seq_record in tqdm(fasta_sequences, desc="Processing Sequences", unit="sequence"):
            # 创建临时FASTA文件，包含当前处理的序列
            temp_fasta_path = "temp_seq.fasta"
            SeqIO.write(seq_record, temp_fasta_path, "fasta")
 
            # 创建BLAST命令行对象
            blastp_cline = NcbiblastpCommandline(
                query=temp_fasta_path, 
                db=db_path, 
                evalue=10.0,  # Expectation value
                outfmt=5, 
                out="blast_results.xml",  # 暂存BLAST结果的XML文件
                matrix="BLOSUM62",  # Weight Matrix
                num_descriptions=5,  # Max Scores
                num_alignments=5,  # Max Alignments（与max_hits保持一致）
                gapopen=11,  # Standard default for BLAST protein search gap opening
                gapextend=1,  # Standard default for BLAST protein search gap extension
                word_size=3  # Word size
            )
            stdout, stderr = blastp_cline()
 
            # 使用with语句来确保文件在解析完毕后被正确关闭
            with open("blast_results.xml") as result_handle:
                blast_records = list(NCBIXML.parse(result_handle))
 
            # 处理当前序列的BLAST结果
            output_file.write(f"Query: {seq_record.id}\n")
            if blast_records:
                hits = sorted(blast_records[0].alignments, key=lambda aln: aln.hsps[0].score, reverse=True)[:max_hits]
 
                for hit in hits:
                    # 使用hit_accession作为输出基因名称
                    output_file.write(f"\tHit: {hit.accession}, Score: {hit.hsps[0].score}\n")
            else:
                output_file.write("\tNo hits found\n")
 
            # 删除临时文件
            os.remove(temp_fasta_path)
            os.remove("blast_results.xml")
 
if __name__ == "__main__":
    fasta_file_path = "要比对的文件.fasta"  # 替换为你的FASTA文件路径
    blast_db_path = "建好的数据库"  # 替换为你的BLAST数据库路径
    output_file_path = "输出的文件"  # 定义输出文件路径
 
    blast_fasta(fasta_file_path, blast_db_path, output_file_path)
