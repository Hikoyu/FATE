#!/usr/bin/env perl

# Copyright (c) 2016-2018 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;
no warnings 'portable';

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "fate.pl";	# ソフトウェアの名前
my $version = "ver.2.0.0";	# ソフトウェアのバージョン
my $note = "FATE is Framework for Annotating Translatable Exons.\n  This software annotates protein-coding regions by a classical homology-based method.";	# ソフトウェアの説明
my $usage = "<required items> [optional items]";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
$command{"search"} = "Search protein-coding regions from genomic DNA sequences under specified conditions";
$command{"filter"} = "Filter already annotated protein-coding regions under specified conditions";
$command{"predict"} = "Predict protein-coding regions from cDNA sequences under specified conditions";
### 編集範囲 終了 ###
my @command_list = sort(keys(%command));

# 指定されたコマンドを確認
my $specified_command = shift(@ARGV) if @command_list and @ARGV;
&exception::error("unknown command: $specified_command") if $specified_command and !grep {$_ eq $specified_command} @command_list;

# 共通オプションを定義
my %option;
### 編集範囲 開始 ###
$option{"w"} = "\tUse 2-byte line feed code (CR+LF) for input files";
$option{"b PATH "} = "Path to blast output file (TSV format)";
$option{"f STR "} = "Output format <bed|gtf> [bed]";
$option{"n STR "} = "Prefix of each locus name [locus]";
$option{"p INT "} = "Number of parallel worker threads <1-> [1]";
$option{"t INT "} = "Output biotype (+1:functional genes, +2:truncated genes, +4:pseudogenes) <1-7> [7]";
$option{"v INT "} = "Maximum number of candidate isoforms <1->";
### 編集範囲 終了 ###

# コマンドごとのオプション定義を取得
&{\&{"${specified_command}::define"}} if $specified_command;
my @option_list = sort(keys(%option));

# ヘルプを表示 (引数未指定時)
&exception::help if !@ARGV and !-p STDIN;

# オプションの入力処理
my %opt;
$_ = join("", @option_list);
$_ =~ s/\s+\S+\s+/:/g;
getopts($_, \%opt);

# 未指定オプションのデフォルト値を入力
foreach (@option_list) {
	$opt{substr($_, 0, 1)} = substr($option{$_}, index($option{$_}, "[") + 1, index($option{$_}, "]") - index($option{$_}, "[") - 1) if $option{$_} =~ /\[.+\]$/ and !defined($opt{substr($_, 0, 1)});
}

### 編集範囲 開始 ###
# 追加のモジュールを宣言
use List::Util;
use IPC::Open2;
use IO::Pipe;
use threads;
use Thread::Queue;

# コドンを定義
my %codon = (
	"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
	"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
	"TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
	"TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
	"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
	"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
	"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
	"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
	"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
	"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
	"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
	"AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
	"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
	"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
	"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
	"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G",
	"TTY" => "F", "TTR" => "L",
	"TCY" => "S", "TCR" => "S", "TCM" => "S", "TCK" => "S", "TCS" => "S", "TCW" => "S", "TCH" => "S", "TCB" => "S", "TCV" => "S", "TCD" => "S", "TCN" => "S", "TC" => "S",
	"TAY" => "Y", "TAR" => "*",
	"TGY" => "C",
	"CTY" => "L", "CTR" => "L", "CTM" => "L", "CTK" => "L", "CTS" => "L", "CTW" => "L", "CTH" => "L", "CTB" => "L", "CTV" => "L", "CTD" => "L", "CTN" => "L", "CT" => "L",
	"CCY" => "P", "CCR" => "P", "CCM" => "P", "CCK" => "P", "CCS" => "P", "CCW" => "P", "CCH" => "P", "CCB" => "P", "CCV" => "P", "CCD" => "P", "CCN" => "P", "CC" => "P",
	"CAY" => "H", "CAR" => "Q",
	"CGY" => "R", "CGR" => "R", "CGM" => "R", "CGK" => "R", "CGS" => "R", "CGW" => "R", "CGH" => "R", "CGB" => "R", "CGV" => "R", "CGD" => "R", "CGN" => "R", "CG" => "R",
	"ATY" => "I",               "ATM" => "I",                             "ATW" => "I", "ATH" => "I",
	"ACY" => "T", "ACR" => "T", "ACM" => "T", "ACK" => "T", "ACS" => "T", "ACW" => "T", "ACH" => "T", "ACB" => "T", "ACV" => "T", "ACD" => "T", "ACN" => "T", "AC" => "T",
	"AAY" => "N", "AAR" => "K",
	"AGY" => "S", "AGR" => "R",
	"GTY" => "V", "GTR" => "V", "GTM" => "V", "GTK" => "V", "GTS" => "V", "GTW" => "V", "GTH" => "V", "GTB" => "V", "GTV" => "V", "GTD" => "V", "GTN" => "V", "GT" => "V",
	"GCY" => "A", "GCR" => "A", "GCM" => "A", "GCK" => "A", "GCS" => "A", "GCW" => "A", "GCH" => "A", "GCB" => "A", "GCV" => "A", "GCD" => "A", "GCN" => "A", "GC" => "A",
	"GAY" => "D", "GAR" => "E",
	"GGY" => "G", "GGR" => "G", "GGM" => "G", "GGK" => "G", "GGS" => "G", "GGW" => "G", "GGH" => "G", "GGB" => "G", "GGV" => "G", "GGD" => "G", "GGN" => "G", "GG" => "G",
	"YTA" => "L", "YTG" => "L", "YTR" => "L",
	"MGA" => "R", "MGG" => "R", "MGR" => "R"
);

# 変換テンプレートを定義
my $bed6_mask_template = "S/A*LLS/A*C/A*cL/a*L/a*";
my $bed12_template = "S/A*LLS/A*C/A*cLLC/A*LL/a*L/a*";

# 処理を追加
### 編集範囲 終了 ###

# メインルーチンを実行
&main;
exit(0);

# メインルーチン
sub main {
	### 編集範囲 開始 ###
	# 指定された共通オプションを確認
	$opt{"f"} = lc($opt{"f"});
	&exception::error("specify bed or gtf: -f $opt{f}") if $opt{"f"} ne "bed" and $opt{"f"} ne "gtf";
	&exception::error("specify INT >= 1: -p $opt{p}") if $opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1;
	&exception::error("specify INT 1-7: -t $opt{t}") if $opt{"t"} ne "1" and $opt{"t"} ne "2" and $opt{"t"} ne "3" and $opt{"t"} ne "4" and $opt{"t"} ne "5" and $opt{"t"} ne "6" and $opt{"t"} ne "7";
	&exception::error("specify INT >= 1: -v $opt{v}") if defined($opt{"v"}) and ($opt{"v"} !~ /^\d+$/ or $opt{"v"} < 1);
	$opt{"v"} = 0 if !defined($opt{"v"});
	
	# 相同性検索の出力ファイルを確認 (-b指定時)
	if ($opt{"b"}) {
		&exception::error("file not found: $opt{b}") if !-f $opt{"b"};
		&exception::error("file unreadable: $opt{b}") if !-r $opt{"b"};
		&exception::error("null file specified: $opt{b}") if !-s $opt{"b"};
	}
	### 編集範囲 終了 ###
	
	# コマンドの実行 (コマンド指定時)
	&{\&{"${specified_command}::body"}} if $specified_command;
	
	### 編集範囲 開始 ###
	# 処理を追加
	### 編集範囲 終了 ###
	return(1);
}

## ここから例外処理のパッケージ ##
package exception;

# ヘルプ表示
sub help {
	print STDERR "$software ";
	print STDERR $specified_command ? $specified_command : $version;
	print STDERR "\n\nFunctions:\n  $note\n\nUsage:\n  $software ";
	if (!$specified_command and @command_list) {
		print STDERR "<command>\n";
		print STDERR "\nCommand:\n";
		foreach (@command_list) {print STDERR "  $_\t$command{$_}\n";}
	}
	else {
		print STDERR "$specified_command " if $specified_command;
		print STDERR "[options] " if @option_list;
		print STDERR "$usage\n";
		print STDERR "\nOptions:\n" if @option_list;
		foreach (@option_list) {print STDERR "  -$_\t$option{$_}\n";}
	}
	exit(0);
}

# エラー表示
sub error {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Error: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	exit(1);
}

# 注意表示
sub caution {
	print STDERR $software;
	print STDERR " $specified_command" if $specified_command;
	print STDERR ": Caution: $_[0]";
	print STDERR ": $_[1] line $." if $_[1];
	print STDERR "\n";
	return(1);
}

### 編集範囲 開始 ###
## ここからsearchコマンドのパッケージ ##
package search;

# コマンドとオプションを定義
sub define {
	$note = "Search protein-coding regions from genomic DNA sequences under specified conditions.";
	$usage = "<genome.fna> <STDIN|in1.fa> [in2.fa ...] [>out.bed|>out.gtf]";
	$option{"s"} = "\tForce GT-AG rule for splice junctions";
	$option{"g STR "} = "Gene prediction engine <exonerate|genewise>";
	$option{"h STR "} = "Homology search engine <blastn|dc-megablast|megablast|tblastn|tblastn-fast> [blastn]";
	$option{"5 INT "} = "Length of 5' flanking region <0-> [300]";
	$option{"3 INT "} = "Length of 3' flanking region <0-> [300]";
	$option{"i INT "} = "Maximum interval length allowed to assemble initial hits [100000]";
	$option{"o INT "} = "Maximum overlap/gap length of query boundries allowed to assemble initial hits [30]";
	$option{"l INT "} = "Minimum length to regard as complete CDS <0-> [0]";
	$option{"c FLOAT "} = "Minimum query coverage to regard as complete CDS <0-1> [0.85]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("unknown engine specified: -g $opt{g}") if $opt{"g"} and $opt{"g"} ne "exonerate" and $opt{"g"} ne "genewise";
	&exception::error("unknown engine specified: -h $opt{h}") if $opt{"h"} ne "blastn" and $opt{"h"} ne "dc-megablast" and $opt{"h"} ne "megablast" and $opt{"h"} ne "tblastn" and $opt{"h"} ne "tblastn-fast";
	&exception::error("specify INT >= 0: -5 $opt{5}") if $opt{"5"} !~ /^\d+$/;
	&exception::error("specify INT >= 0: -3 $opt{3}") if $opt{"3"} !~ /^\d+$/;
	&exception::error("specify INT >= 0: -i $opt{i}") if $opt{"i"} !~ /^\d+$/;
	&exception::error("specify INT >= 0: -o $opt{o}") if $opt{"o"} !~ /^\d+$/;
	&exception::error("specify INT >= 0: -l $opt{l}") if $opt{"l"} !~ /^\d+$/;
	&exception::error("specify FLOAT 0-1: -c $opt{c}") if $opt{"c"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"c"} > 1;
	&exception::caution("-s ignored under -g unspecified") if $opt{"s"} and !$opt{"g"};
	
	# 依存関係を確認
	my $software_list = &common::check_dependencies(["makeblastdb", "blastn", "tblastn", $opt{"g"}]);
	&exception::error("makeblastdb not installed") if !$software_list->{"makeblastdb"};
	&exception::error("blastn not installed") if !$software_list->{"blastn"};
	&exception::error("tblastn not installed") if !$software_list->{"tblastn"};
	&exception::error("$opt{g} not installed") if $opt{"g"} and !$software_list->{$opt{"g"}};
	
	# オプションの対応を確認
	&exception::error("$opt{h} unavailable on your blast+: -h $opt{h}") if $opt{"h"} eq "tblastn-fast" and $software_list->{"tblastn"} =~ /blast 2\.2\.[1-2]\d/;
	
	# ゲノム配列のfastaファイル名を取得
	my $genome_file = shift(@ARGV);
	
	# 入力ファイルを確認
	&exception::error("input file not specified") if !@ARGV and !-p STDIN;
	foreach (@ARGV) {
		&exception::error("file not found: $_") if !-f $_;
		&exception::error("file unreadable: $_") if !-r $_;
		&exception::error("null file specified: $_") if !-s $_;
	}
	
	# ゲノム配列のfastaファイルを検索
	&common::find_db($genome_file);
	
	# ゲノム配列のfastaインデックスを取得
	my $genome_faidx = &common::read_fasta($genome_file, $opt{"w"});
	
	# ゲノム配列のblastデータベースを確認 (-b未指定時)
	&common::check_blastdb($genome_file) if !$opt{"b"};
	
	# デフォルトの処理を定義
	my %homology_search = (
		"blastn"       => "blastn  -num_threads $opt{p} -db $genome_file -outfmt 6 -soft_masking true -task blastn",
		"dc-megablast" => "blastn  -num_threads $opt{p} -db $genome_file -outfmt 6 -soft_masking true -task dc-megablast",
		"megablast"    => "blastn  -num_threads $opt{p} -db $genome_file -outfmt 6 -soft_masking true -task megablast",
		"tblastn"      => "tblastn -num_threads $opt{p} -db $genome_file -outfmt 6 -soft_masking true",
		"tblastn-fast" => "tblastn -num_threads $opt{p} -db $genome_file -outfmt 6 -soft_masking true -task tblastn-fast"
	);
	my %gene_prediction = (
		"exonerate" => "exonerate -m protein2genome -V 0 --showalignment F --showvulgar F --showtargetgff T -r F",
		"genewise"  => "genewise -divide '#' -pseudo -sum -gff -silent -quiet"
	);
	
	# オプションの処理を追加
	foreach (keys(%homology_search)) {$homology_search{$_} .= `which tee` ? " 2>/dev/null | tee fate_search_$opt{h}.out" : "";}
	$gene_prediction{"exonerate"} .= $opt{"s"} ? " --forcegtag T" : "";
	$gene_prediction{"genewise"} .= $opt{"s"} ? "" : " -nosplice_gtag";
	
	# 個々の配列データを保存しておくディレクトリを作成 (-g指定時)
	mkdir("prot") or &exception::error("failed to make directory: prot") if $opt{"g"} and !-d "prot";
	mkdir("nucl") or &exception::error("failed to make directory: nucl") if $opt{"g"} and !-d "nucl";
	
	# 相同性検索の出力ファイルを開く (-b指定時)
	open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}") if $opt{"b"};
	
	# 相同性検索を実行 (-b未指定時)
	main::open2(*SEARCH_OUT, *SEARCH_IN, "$homology_search{$opt{h}} 2>/dev/null") or &exception::error("failed to execute homology search") if !$opt{"b"};
	
	# プロセス間通信のパイプを作成
	my $report = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	&exception::error("failed to fork process") if !defined($pid);
	
	## ここから子プロセスの処理 ##
	if (!$pid) {
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# 変数を宣言
		my $query_title = undef;
		my $query_seq = "";
		
		# パイプを開く
		$report->writer;
		
		# 入力の改行コードを一時的に変更 (-w指定時)
		local $/ = "\r\n" if $opt{"w"};
		
		# クエリー配列を読み込みながら処理
		while (my $line = <>) {
			# 改行コードを除去
			chomp($line);
			
			# 配列行の処理
			$query_seq .= $line =~ /^>/ ? "" : uc($line);
			
			# ID行の処理
			if ($line =~ /^>/ or eof) {
				# 配列データを処理
				if (defined($query_title)) {
					# 親プロセスにクエリー配列長を送信
					print $report length($query_seq), "\n";
					
					# クエリー配列をファイルに出力 (-g指定時)
					if ($opt{"g"}) {
						# 変数を宣言
						my $aa_seq = $opt{"h"} =~ /^tblastn/ ? $query_seq : &common::translate($query_seq, 0);
						
						# 配列中に含まれるアスタリスクを除去
						$aa_seq =~ s/\*//g;
						
						# 出力ファイルを作成
						open(PROT, ">", "prot/$query_title.fa") or &exception::error("failed to make file: prot/$query_title.fa");
						
						# クエリー配列をfasta形式でファイルに出力
						print PROT ">$query_title\n$aa_seq\n";
						
						# 出力ファイルを閉じる
						close(PROT);
					}
					
					# クエリー配列を相同性検索に入力 (-b未指定時)
					print SEARCH_IN ">$query_title\n$query_seq\n" if !$opt{"b"};
				}
				
				# ID行の最初の空白文字の前までを配列名として登録
				($query_title) = split(/\s/, substr($line, 1));
				
				# 配列をリセット
				$query_seq = "";
			}
		}
		
		# パイプを閉じる
		close($report);
		
		# 相同性検索への入力を閉じる (-b未指定時)
		close(SEARCH_IN) if !$opt{"b"};
		
		# プロセスを終了
		exit(0);
	}
	## ここまで子プロセスの処理 ##
	
	# 相同性検索への入力を閉じる (-b未指定時)
	close(SEARCH_IN) if !$opt{"b"};
	
	# 変数を宣言
	my %assembled_hits = ();
	my %query_len = ();
	my @blast_hits = ();
	my $last_query = undef;
	my $last_subject = "";
	my $num_assembly = 0;
	
	# パイプを開く
	$report->reader;
	
	# 相同性検索の出力を読み込みながら処理
	print STDERR $opt{"b"} ? "Loading homology search results..." : "Running homology search...";
	while (my $line = <SEARCH_OUT>) {
		# コメント行を除外
		next if $line =~ /^["#"]/;
		
		# 改行コードを除去
		chomp($line);
		
		# タブ文字でデータを分割
		my @col = split(/\t/, $line);
		
		# クエリー名またはサブジェクト名が変わった場合
		if (!defined($last_query) or $col[0] ne $last_query or $col[1] ne $last_subject) {
			# ここまでのヒットをアセンブル
			$assembled_hits{$last_query}{$last_subject} = &assemble(\@blast_hits, $opt{"i"}, $opt{"o"}) and $num_assembly += scalar(map {keys(%{$_->{"assemble"}})} @{$assembled_hits{$last_query}{$last_subject}}) if defined($last_query);
			
			# クエリー名が変わった場合は子プロセスからクエリー配列長を受信し、クエリー名を更新
			$query_len{$col[0]} .= <$report> and $last_query = $col[0] if !defined($last_query) or $col[0] ne $last_query;
			
			# サブジェクト名を更新
			$last_subject = $col[1];
			
			# ヒットリストをリセット
			@blast_hits = ();
		}
		
		# 検索結果をヒットリストに追加
		push(@blast_hits, $col[8] < $col[9] ? {"query_start" => $col[6] - 1, "query_end" => $col[7], "locus_start" => $col[8] - 1, "locus_end" => $col[9], "strand" => 1} : {"query_start" => $col[7], "query_end" => $col[6] - 1, "locus_start" => $col[9] - 1, "locus_end" => $col[8], "strand" => -1});
		$blast_hits[-1]->{"assemble"} = {};
		$blast_hits[-1]->{"score"} = $col[11];
		$blast_hits[-1]->{"num_connection"} = 0;
	}
	
	# 残りのヒットをアセンブル
	$assembled_hits{$last_query}{$last_subject} = &assemble(\@blast_hits, $opt{"i"}, $opt{"o"}) and $num_assembly += scalar(map {keys(%{$_->{"assemble"}})} @{$assembled_hits{$last_query}{$last_subject}});
	
	# パイプを閉じる
	close($report);
	
	# 相同性検索の出力を閉じる
	close(SEARCH_OUT);
	
	# 子プロセスが終了するまで待機
	waitpid($pid, 0);
	
	# 子プロセスが異常終了した場合
	&exception::error("process abnormally exited") if $?;
	print STDERR "completed\n";
	
	# 相同性検索でヒットが得られなかった場合
	&exception::error("no hits found from $opt{h} search") if !defined($last_query);
	
	# 変数を宣言
	my %loci = ();
	
	# 遺伝子構造予測を実施する場合 (-g指定時)
	if ($opt{"g"}) {
		# 入出力のキューを作成
		my $input = Thread::Queue->new;
		my $output = Thread::Queue->new;
		
		# 入力キューにアセンブルデータをバイナリ形式で追加
		print STDERR "Running gene prediction...";
		foreach my $query (keys(%assembled_hits)) {
			foreach my $subject (keys(%{$assembled_hits{$query}})) {
				foreach my $locus (@{$assembled_hits{$query}{$subject}}) {
					foreach my $assemble (values(%{$locus->{"assemble"}})) {
						$input->enqueue(pack($bed6_mask_template, $subject, $locus->{"locus_start"}, $assemble->{"locus_destination"}, $query, $opt{"h"} =~ /^tblastn/ ? $query_len{$query} * 3 : $query_len{$query}, $locus->{"strand"}, pack("L*", @{$assemble->{"mask_block_size"}}), pack("L*", @{$assemble->{"mask_block_start"}})));
					}
				}
			}
		}
		print STDERR "\rRunning gene prediction...0%";
		
		# 変数を宣言
		my $num_data = $input->pending;
		
		# 指定したワーカースレッド数で並列処理
		for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
			## ここからワーカースレッドの処理 ##
			threads::async {
				# ワーカースレッドのみを終了可能に変更
				threads->set_thread_exit_only(1);
				
				# ゲノム配列のfastaファイルを開く
				open(GENOME, "<", $genome_file) or &exception::error("failed to open file: $genome_file");
				
				# 入力キューからデータをバイナリ形式で取得して処理
				while (defined(my $dat = $input->dequeue_nb)) {
					# データを変換
					my ($subject, $locus_start, $locus_end, $query, $query_len, $strand, $mask_block_size, $mask_block_start) = unpack($bed6_mask_template, $dat);
					
					# 変数を宣言
					my $faidx = &common::faidx_decode($genome_faidx->{$subject});
					
					# マスクブロックサイズとマスクブロック開始点リストを取得
					my @mask_block_size = unpack("L*", $mask_block_size);
					my @mask_block_start = unpack("L*", $mask_block_start);
					
					# フランキング配列長に合わせて開始点と終了点を修正
					my $modified_locus_start = $strand > 0 ? $locus_start - $opt{"5"} : $locus_start - $opt{"3"};
					my $modified_locus_end = $strand > 0 ? $locus_end + $opt{"3"} : $locus_end + $opt{"5"};
					$modified_locus_start = 0 if $modified_locus_start < 0;
					$modified_locus_end = $faidx->{"seq_length"} if $modified_locus_end > $faidx->{"seq_length"};
					
					# 取り出す領域の位置情報を算出
					my $start_point = int($modified_locus_start / $faidx->{"row_width"}) * $faidx->{"row_bytes"} + $modified_locus_start % $faidx->{"row_width"};
					my $end_point = int($modified_locus_end / $faidx->{"row_width"}) * $faidx->{"row_bytes"} + $modified_locus_end % $faidx->{"row_width"};
					
					# 領域の塩基配列を取得
					seek(GENOME, $faidx->{"seq_start"} + $start_point, 0);
					read(GENOME, my $locus_seq, $end_point - $start_point);
					
					# 配列中に含まれる改行コードを除去
					$locus_seq =~ s/\n|\r//g;
					
					# 配列を大文字に変換
					$locus_seq = uc($locus_seq);
					
					# 配列をマスク
					for (my $i = 0;$i < @mask_block_size;$i++) {substr($locus_seq, $mask_block_start[$i] - $modified_locus_start, $mask_block_size[$i]) = "N" x $mask_block_size[$i];}
					
					# 相補鎖に変換
					&common::complementary($locus_seq) if $strand < 0;
					
					# 領域開始点を表示形式に合わせる
					$modified_locus_start++;
					
					# 領域名を定義
					my $locus_name = "$subject:$modified_locus_start-$modified_locus_end";
					$locus_name .= $strand > 0 ? "(+)" : "(-)";
					
					# 領域の配列ファイルを作成
					open(NUCL, ">", "nucl/locus$thread_id.fa") or &exception::error("failed to make file: nucl/locus$thread_id.fa");
					
					# 取得した配列をファイルに出力
					print NUCL ">$locus_name\n$locus_seq\n";
					
					# 領域の配列ファイルを閉じる
					close(NUCL);
					
					# 変数を宣言
					my @genes = ();
					my $query_file = "";
					my $target_file = "";
					my $query_start = 0;
					my $query_end = 0;
					my $summary_flag = 0;
					
					# 遺伝子構造予測の引数を定義
					$query_file = "-q prot/$query.fa" if $opt{"g"} eq "exonerate";
					$query_file = "prot/$query.fa" if $opt{"g"} eq "genewise";
					$target_file = "-t nucl/locus$thread_id.fa" if $opt{"g"} eq "exonerate";
					$target_file = "nucl/locus$thread_id.fa" if $opt{"g"} eq "genewise";
					$query_file =~ s/\|/\\\|/g;
					
					# 遺伝子構造予測を実行
					open(PREDICT, "-|", "$gene_prediction{$opt{g}} $query_file $target_file 2>/dev/null") or &exception::error("failed to execute gene prediction: $query vs $subject");
					
					# 領域開始点をbed形式に戻す
					$modified_locus_start--;
					
					# 遺伝子構造予測結果を1行ずつ読み込んで処理
					while (my $line = <PREDICT>) {
						# コメント行を除外
						++$summary_flag and next if $line =~ /^["#"]/;
						
						# 空白文字でデータを分割
						my @col = split(/\s+/, $line, 9);
						
						# 見出し行を除外
						next if $col[0] eq "Bits" and $col[1] eq "Query" and $col[4] eq "Target";
						
						# summary行の処理 (-g genewise指定時)
						($query_start, $query_end) = (($col[2] - 1) * 3, $col[3] * 3) and next if $opt{"g"} eq "genewise" and !($summary_flag & 0x01);
						
						# 開始点をbed形式に合わせる
						$col[3]--;
						
						# match行またはgene行の処理
						push(@genes, [$subject, $col[3], $col[4], $query, $col[5], $strand, 0, 0, ".", 0, [], [], $query_start, $query_end]) if $col[2] eq "match" or $col[2] eq "gene";
						
						# cds行の処理
						++$genes[-1]->[9] and push(@{$genes[-1]->[10]}, $col[4] - $col[3]) and push(@{$genes[-1]->[11]}, $col[3]) if $col[2] eq "cds";
						
						# similarity行の処理 (-g exonerate指定時)
						my @query_block = map {[split(/\s/)]} grep {/^Align/} split(/ ; /, $col[8]) if $opt{"g"} eq "exonerate" and $col[2] eq "similarity";
						$genes[-1]->[12] = List::Util::min(map {($_->[2] - 1) * 3} @query_block) if @query_block;
						$genes[-1]->[13] = List::Util::max(map {($_->[2] - 1) * 3 + $_->[3]} @query_block) if @query_block;
					}
					
					# 遺伝子構造予測を終了
					close(PREDICT);
					
					# 各遺伝子について翻訳領域を推定
					foreach my $gene (@genes) {
						# 変数を宣言
						my $upstream_truncation = 0;
						my $downstream_truncation = 0;
						my $completeness = 0;
						
						# Nを含まない上流配列を取得し、その長さが指定値より短い場合は5'側truncateとする
						my $upstream_seq = substr(substr($locus_seq, 0, $gene->[1]), -$opt{"5"});
						my $upstream_offset = $strand > 0 ? $locus_start - $modified_locus_start - $gene->[1] : $modified_locus_end - $gene->[1] - $locus_end;
						$upstream_offset = 0 if $upstream_offset < 0;
						$upstream_seq =~ s/^.*N//;
						$upstream_truncation = 1 if length($upstream_seq) < $opt{"5"} - $upstream_offset;
						
						# Nを含まない下流配列を取得し、その長さが指定値より短い場合は3'側truncateとする
						my $downstream_seq = substr($locus_seq, $gene->[2], $opt{"3"});
						my $downstream_offset = $strand > 0 ? $modified_locus_start + $gene->[2] - $locus_end : $locus_start - $modified_locus_end + $gene->[2];
						$downstream_offset = 0 if $downstream_offset < 0;
						$downstream_seq =~ s/N.*$//;
						$downstream_truncation = 1 if length($downstream_seq) < $opt{"3"} - $downstream_offset;
						
						# 5'側クエリー被覆率が指定値以上の場合
						if ($query_len - $gene->[12] >= $query_len * $opt{"c"}) {
							# 上流フランキング領域のアミノ酸配列を取得
							my $upstream_aa = &common::translate($upstream_seq, length($upstream_seq) % 3);
							
							# 先頭エキソンのアミノ酸配列を取得
							my $first_exon_aa = &common::translate(substr($locus_seq, $gene->[11]->[0], $gene->[10]->[0]), 0);
							
							# 開始コドンを探索
							my $start_pos = rindex($upstream_aa, "M");
							my $stop_pos = rindex($upstream_aa, "*");
							$start_pos = $stop_pos < $start_pos ? $start_pos - length($upstream_aa) : index($first_exon_aa, "M") + 1;
							$start_pos -= $start_pos > 0;
							$start_pos *= 3;
							
							# 開始コドンを考慮しても5'側クエリー被覆率が指定値以上の場合は先頭ブロックを修正
							($gene->[10]->[0], $gene->[11]->[0]) = ($gene->[10]->[0] - $start_pos, $gene->[11]->[0] + $start_pos) if $query_len - $gene->[12] - $start_pos >= $query_len * $opt{"c"};
							
							# 5'端が完全であることを認定
							$upstream_truncation = 0;
							$completeness++;
						}
						
						# 3'側クエリー被覆率が指定値以上の場合
						if ($gene->[13] >= $query_len * $opt{"c"}) {
							# 下流フランキング領域のアミノ酸配列を取得
							my $downstream_aa = &common::translate($downstream_seq, 0);
							
							# 末尾エキソンのアミノ酸配列を取得
							my $last_exon_aa = &common::translate(substr($locus_seq, $gene->[11]->[-1], $gene->[10]->[-1]), $gene->[10]->[-1] % 3);
							
							# 終止コドンを探索
							my $terminal_pos = rindex($last_exon_aa, "*");
							$terminal_pos = $terminal_pos < 0 ? index($downstream_aa, "*") : $terminal_pos - length($last_exon_aa);
							$terminal_pos++;
							$terminal_pos *= 3;
							
							# 終止コドンを考慮しても3'側クエリー被覆率が指定値以上の場合は末尾ブロックを修正
							$gene->[10]->[-1] += $terminal_pos if $gene->[13] + $terminal_pos >= $query_len * $opt{"c"};
							
							# 3'端が完全であることを認定
							$downstream_truncation = 0;
							$completeness++;
						}
						
						# 領域開始点・終了点を修正
						($gene->[1], $gene->[2]) = ($gene->[11]->[0], $gene->[10]->[-1] + $gene->[11]->[-1]);
						
						# ブロックの基準点を修正
						my $basal_pos = $gene->[11]->[0];
						
						# ブロックの相対位置を修正
						$gene->[11] = [map {$_ - $basal_pos} @{$gene->[11]}];
						
						# 翻訳領域のアミノ酸配列を取得
						my $cds = "";
						for (my $i = 0;$i < $gene->[9];$i++) {$cds .= substr($locus_seq, $gene->[1] + $gene->[11]->[$i], $gene->[10]->[$i]);}
						my $aa = &common::translate($cds, 0);
						
						# フレームシフトが存在する場合 (偽遺伝子)
						if (List::Util::sum(@{$gene->[10]}) % 3 > 0) {$gene->[8] = "red";}
						
						# 終止コドンが末尾以外に存在する場合 (偽遺伝子)
						elsif (index($aa, "*") >= 0 and index($aa, "*") < length($aa) - 1) {$gene->[8] = "red";}
						
						# 両末端が完全とみなされない場合 (分断遺伝子または偽遺伝子)
						elsif ($completeness < 2) {$gene->[8] = $upstream_truncation | $downstream_truncation ? "yellow" : "red";}
						
						# 配列長が指定値未満の場合 (偽遺伝子)
						elsif (List::Util::sum(@{$gene->[10]}) < $opt{"l"}) {$gene->[8] = "red";}
						
						# 上記に該当しない場合 (機能遺伝子)
						else {$gene->[8] = "blue";}
						
						# ゲノム配列の座標でデータを修正
						($gene->[1], $gene->[2]) = $strand > 0 ? ($modified_locus_start + $gene->[1], $modified_locus_start + $gene->[2]) : ($modified_locus_end - $gene->[2], $modified_locus_end - $gene->[1]);
						($gene->[6], $gene->[7]) = ($gene->[1], $gene->[2]);
						if ($strand < 0) {
							$gene->[10] = [reverse(@{$gene->[10]})];
							$gene->[11] = [reverse(@{$gene->[11]})];
							for (my $i = 0;$i < @{$gene->[11]};$i++) {$gene->[11]->[$i] = $gene->[2] - $gene->[1] - $gene->[10]->[$i] - $gene->[11]->[$i];}
						}
						
						# データをバイナリ形式に変換
						$gene = pack($bed12_template, @{$gene}[0..9], pack("L*", @{$gene->[10]}), pack("L*", @{$gene->[11]}));
					}
					
					# 出力キューにデータをバイナリ形式で追加
					$output->enqueue(pack("S/A*L(L/a*)*", $subject, scalar(@genes), @genes));
					
					# 完了データ数を取得
					my $fin_data = $output->pending;
					
					# 経過を表示
					print STDERR "\rRunning gene prediction...", int($fin_data / $num_data * 100), "%" if int($fin_data / $num_data * 100) > int(($fin_data - 1) / $num_data * 100);
				}
				
				# ゲノム配列のfastaファイルを閉じる
				close(GENOME);
				
				# スレッドを終了
				return(1);
			};
			## ここまでワーカースレッドの処理 ##
		}
		
		# 変数を宣言
		my $fin_threads = 0;
		
		# 各ワーカースレッドが終了するまで待機
		foreach (threads->list) {$_->join and $fin_threads++;}
		
		# ワーカースレッドが異常終了した場合
		&exception::error("thread abnormally exited") if $fin_threads < $opt{"p"};
		
		# 出力キューからデータをバイナリ形式で取得して処理
		while (my $dat = $output->dequeue_nb) {
			# データを変換
			my ($subject, @locus) = unpack("S/A*L/(L/a*)*", $dat);
			
			# ハッシュにデータをバイナリ形式で登録
			push(@{$loci{$subject}}, @locus) if @locus;
		}
		print STDERR "\rRunning gene prediction...completed\n";
	}
	
	# 遺伝子構造予測を実施しない場合 (-g未指定時)
	else {
		# ハッシュにアセンブルデータをバイナリ形式で登録
		print STDERR "Loading assembled loci...";
		foreach my $query (keys(%assembled_hits)) {
			foreach my $subject (keys(%{$assembled_hits{$query}})) {
				foreach my $locus (@{$assembled_hits{$query}{$subject}}) {
					foreach my $assemble (values(%{$locus->{"assemble"}})) {
						push(@{$loci{$subject}}, pack($bed12_template, $subject, $locus->{"locus_start"}, $assemble->{"locus_destination"}, $query, $assemble->{"total_score"}, $locus->{"strand"}, $locus->{"locus_start"}, $assemble->{"locus_destination"}, ".", scalar(@{$assemble->{"block_size"}}), pack("L*", @{$assemble->{"block_size"}}), pack("L*", map {$_ - $locus->{"locus_start"}} @{$assemble->{"block_start"}})));
					}
				}
			}
		}
		print STDERR "completed\n";
	}
	
	# 入出力のキューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入力キューにアセンブルデータをバイナリ形式で追加
	print STDERR "Running isoform definition...";
	foreach my $subject (keys(%loci)) {$input->enqueue(pack("S/A*L(L/a*)*", $subject, scalar(@{$loci{$subject}}), @{$loci{$subject}}));}
	print STDERR "\rRunning isoform definition...0%";
	
	# 変数を宣言
	my $num_data = $input->pending;
	
	# 指定したワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		threads::async {
			# ワーカースレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューからデータをバイナリ形式で取得して処理
			while (defined(my $dat = $input->dequeue_nb)) {
				# データを変換
				my ($subject, @locus) = unpack("S/A*L/(L/a*)*", $dat);
				
				# アイソフォームを分離
				my $isoforms = &common::define_isoforms(\@locus, $opt{"v"}, $opt{"t"});
				
				# アイソフォームを分離し、出力キューにデータをバイナリ形式で追加
				$output->enqueue(pack("S/A*L(L/a*)*", $subject, scalar(@{$isoforms}), @{$isoforms}));
				
				# 完了データ数を取得
				my $fin_data = $output->pending;
				
				# 経過を表示
				print STDERR "\rRunning isoform definition...", int($fin_data / $num_data * 100), "%" if int($fin_data / $num_data * 100) > int(($fin_data - 1) / $num_data * 100);
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 変数を宣言
	my $fin_threads = 0;
	
	# 各ワーカースレッドが終了するまで待機
	foreach (threads->list) {$_->join and $fin_threads++;}
	
	# ワーカースレッドが異常終了した場合
	&exception::error("thread abnormally exited") if $fin_threads < $opt{"p"};
	
	# 変数を宣言
	my $gene_id = 0;
	
	# 出力キューからデータをバイナリ形式で取得して処理
	while (my $dat = $output->dequeue_nb) {
		# データを変換
		my ($subject, @locus) = unpack("S/A*L/(L/a*)*", $dat);
		
		# 領域データリストにデータをバイナリ形式で追加
		$loci{$subject} = \@locus;
	}
	
	# 結果を出力
	foreach my $subject (sort {&common::faidx_decode($genome_faidx->{$a})->{"id_order"} <=> &common::faidx_decode($genome_faidx->{$b})->{"id_order"}} keys(%loci)) {&common::output($loci{$subject}, $gene_id, $opt{"n"}, $opt{"f"});}
	print STDERR "\rRunning isoform definition...completed\n";
	return(1);
}

# 相同性検索ヒットのアセンブル search::assemble(相同性検索ヒットリストリファレンス, 最大ヒット間距離, 最大許容クエリー境界オーバーラップ/ギャップ)
sub assemble {
	# 変数を宣言
	my $assemble_id = 0;
	
	# ヒットを領域開始点と領域終了点の順で並べ替え
	my @sorted_hits = sort {$a->{"locus_start"} <=> $b->{"locus_start"} or $a->{"locus_end"} <=> $b->{"locus_end"}} @{$_[0]};
	
	# 後方のヒットから判定
	for (my $i = -1;$i >= -@sorted_hits;$i--) {
		# 連結するヒットを決定
		for (my $k = $i + 1;$k < 0;$k++) {
			# 指定値を超える距離のヒットに到達した時点で終了
			last if $sorted_hits[$k]->{"locus_start"} - $sorted_hits[$i]->{"locus_end"} > $_[1];
			
			# 方向が異なるヒットを除外
			next if !($sorted_hits[$k]->{"strand"} + $sorted_hits[$i]->{"strand"});
			
			# 指定値を超えるクエリーオーバーラップのヒットを除外
			next if abs($sorted_hits[$k]->{"query_start"} - $sorted_hits[$i]->{"query_end"}) > $_[2];
			
			# クエリー開始点が前進しないヒットを除外
			next if $sorted_hits[$k]->{"query_start"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_start"} * $sorted_hits[$i]->{"strand"};
			
			# クエリー終了点が前進しないヒットを除外
			next if $sorted_hits[$k]->{"query_end"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_end"} * $sorted_hits[$i]->{"strand"};
			
			# マスクするヒットを決定
			my @mask_block = grep {$_->{"strand"} + $sorted_hits[$i]->{"strand"}} @sorted_hits[$i + 1..$k - 1];
			my @mask_block_size = map {$_->{"locus_end"} - $_->{"locus_start"}} @mask_block;
			my @mask_block_start = map {$_->{"locus_start"}} @mask_block;
			
			# 各アセンブルIDについて処理
			foreach (keys(%{$sorted_hits[$k]->{"assemble"}})) {
				# 総スコアが更新されない同一IDのアセンブルを除外
				next if exists($sorted_hits[$i]->{"assemble"}->{$_}) and $sorted_hits[$i]->{"assemble"}->{$_}->{"total_score"} >= $sorted_hits[$i]->{"score"} + $sorted_hits[$k]->{"assemble"}->{$_}->{"total_score"};
				
				# データを更新
				$sorted_hits[$i]->{"assemble"}->{$_} = {
					"total_score" => $sorted_hits[$i]->{"score"} + $sorted_hits[$k]->{"assemble"}->{$_}->{"total_score"},
					"locus_destination" => $sorted_hits[$k]->{"assemble"}->{$_}->{"locus_destination"},
					"block_size" => [$sorted_hits[$i]->{"locus_end"} - $sorted_hits[$i]->{"locus_start"}, @{$sorted_hits[$k]->{"assemble"}->{$_}->{"block_size"}}],
					"block_start" => [$sorted_hits[$i]->{"locus_start"}, @{$sorted_hits[$k]->{"assemble"}->{$_}->{"block_start"}}],
					"mask_block_size" => [@mask_block_size, @{$sorted_hits[$k]->{"assemble"}->{$_}->{"mask_block_size"}}],
					"mask_block_start" => [@mask_block_start, @{$sorted_hits[$k]->{"assemble"}->{$_}->{"mask_block_start"}}]
				};
			}
			$sorted_hits[$k]->{"num_connection"}++;
		}
		
		# 末尾ヒットの場合はデータを新規登録し、アセンブルIDを更新
		$sorted_hits[$i]->{"assemble"}->{$assemble_id} = {
			"total_score" => $sorted_hits[$i]->{"score"},
			"locus_destination" => $sorted_hits[$i]->{"locus_end"},
			"block_size" => [$sorted_hits[$i]->{"locus_end"} - $sorted_hits[$i]->{"locus_start"}],
			"block_start" => [$sorted_hits[$i]->{"locus_start"}],
			"mask_block_size" => [],
			"mask_block_start" => []
		} and $assemble_id++ if !%{$sorted_hits[$i]->{"assemble"}};
	}
	
	# アセンブルデータリストのリファレンスを返す
	return([grep {!$_->{"num_connection"}} @sorted_hits]);
}

## ここからfilterコマンドのパッケージ ##
package filter;

# コマンドとオプションを定義
sub define {
	$note = "Filter already annotated protein-coding regions under specified conditions.";
	$usage = "<genome.fna> <STDIN|in1.bed> [in2.bed ...] [>out.bed|>out.gtf]";
	$option{"d PATH "} = "Path to gene/protein database file (FASTA format)";
	$option{"h STR "} = "Homology search engine <blastn|dc-megablast|megablast|blastx|blastx-fast> [blastn]";
	$option{"k STR "} = "Keywords for filtering (AND[&], OR[;], BUT[!])";
	$option{"r INT "} = "Cutoff rank of hits <1->";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("options incompatible: -b and -d") if $opt{"b"} and $opt{"d"};
	&exception::error("unknown engine specified: -h $opt{h}") if $opt{"h"} ne "blastn" and $opt{"h"} ne "dc-megablast" and $opt{"h"} ne "megablast" and $opt{"h"} ne "blastx" and $opt{"h"} ne "blastx-fast";
	&exception::error("specify INT >= 1: -r $opt{r}") if defined($opt{"r"}) and ($opt{"r"} !~ /^\d+$/ or $opt{"r"} < 1);
	
	# 依存関係を確認
	my $software_list = &common::check_dependencies(["makeblastdb", "blastn", "blastx"]);
	&exception::error("makeblastdb not installed") if !$software_list->{"makeblastdb"};
	&exception::error("blastn not installed") if !$software_list->{"blastn"};
	&exception::error("blastx not installed") if !$software_list->{"blastx"};
	
	# オプションの対応を確認
	&exception::error("$opt{h} unavailable on your blast+: -h $opt{h}") if $opt{"h"} eq "blastx-fast" and $software_list->{"blastx"} =~ /blast 2\.2\.[1-2]\d/;
	
	# ゲノム配列のfastaファイル名を取得
	my $genome_file = shift(@ARGV);
	
	# 入力ファイルを確認
	&exception::error("input file not specified") if !@ARGV and !-p STDIN;
	foreach (@ARGV) {
		&exception::error("file not found: $_") if !-f $_;
		&exception::error("file unreadable: $_") if !-r $_;
		&exception::error("null file specified: $_") if !-s $_;
	}
	
	# ゲノム配列のfastaファイルを検索
	&common::find_db($genome_file);
	
	# ゲノム配列のfastaインデックスを取得
	my $genome_faidx = &common::read_fasta($genome_file, $opt{"w"});
	
	# 参照配列のfastaファイルを検索し、参照配列のblastデータベースを確認 (-d指定時)
	&common::find_db($opt{"d"}) and &common::check_blastdb($opt{"d"}, $opt{"h"} =~ /^blastx/) if $opt{"d"};
	
	# デフォルトの処理を定義 (-d指定時)
	my %homology_search = (
		"blastn"       => "blastn -num_threads $opt{p} -db $opt{d} -outfmt '6 std salltitles' -soft_masking true -strand plus -task blastn",
		"dc-megablast" => "blastn -num_threads $opt{p} -db $opt{d} -outfmt '6 std salltitles' -soft_masking true -strand plus -task dc-megablast",
		"megablast"    => "blastn -num_threads $opt{p} -db $opt{d} -outfmt '6 std salltitles' -soft_masking true -strand plus -task megablast",
		"blastx"       => "blastx -num_threads $opt{p} -db $opt{d} -outfmt '6 std salltitles' -soft_masking true -strand plus",
		"blastx-fast"  => "blastx -num_threads $opt{p} -db $opt{d} -outfmt '6 std salltitles' -soft_masking true -strand plus -task blastx-fast"
	) if $opt{"d"};
	
	# オプションの処理を追加
	foreach (keys(%homology_search)) {$homology_search{$_} .= `which tee` ? " 2>/dev/null | tee fate_filter_$opt{h}.out" : "";}
	
	# 相同性検索の出力ファイルを開く (-b指定時)
	open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}") if $opt{"b"};
	
	# 相同性検索を実行 (-d指定時)
	main::open2(*SEARCH_OUT, *SEARCH_IN, "$homology_search{$opt{h}} 2>/dev/null") or &exception::error("failed to execute homology search") if $opt{"d"};
	
	# プロセス間通信のパイプを作成
	my $report = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	&exception::error("failed to fork process") if !defined($pid);
	
	## ここから子プロセスの処理 ##
	if (!$pid) {
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# パイプを開く
		$report->writer;
		
		# ゲノム配列のfastaファイルを開く
		open(GENOME, "<", $genome_file) or &exception::error("failed to open file: $genome_file");
		
		# 入力の改行コードを一時的に変更 (-w指定時)
		local $/ = "\r\n" if $opt{"w"};
		
		# bedデータを読み込みながら処理
		while (my $line = <>) {
			# コメント行を除外
			next if $line =~ /^["#"]/;
			
			# 改行コードを除去
			chomp($line);
			
			# タブ文字でデータを分割
			my @col = split(/\t/, $line);
			
			# 変数を宣言
			my $faidx = &common::faidx_decode($genome_faidx->{$col[0]});
			
			# bedデータを確認
			&common::check_bed(\@col, "$opt{n}$.");
			
			# 親プロセスにbedデータを送信 (-d未指定時)
			print $report join("\t", @col), "\n";
			
			# 相同性検索を実行しない場合
			next if !$opt{"d"};
			
			# 取り出す領域の位置情報を算出
			my $start_point = int($col[1] / $faidx->{"row_width"}) * $faidx->{"row_bytes"} + $col[1] % $faidx->{"row_width"};
			my $end_point = int($col[2] / $faidx->{"row_width"}) * $faidx->{"row_bytes"} + $col[2] % $faidx->{"row_width"};
			
			# 領域の塩基配列を取得
			seek(GENOME, $faidx->{"seq_start"} + $start_point, 0);
			read(GENOME, my $locus_seq, $end_point - $start_point);
			
			# 配列中に含まれる改行コードを除去
			$locus_seq =~ s/\n|\r//g;
			
			# 変数を宣言
			my $query_seq = "";
			
			# データを変換
			my @block_size = split(/,/, $col[10]);
			my @block_start = split(/,/, $col[11]);
			
			# ブロックを連結
			for (my $i = 0;$i < $col[9];$i++) {$query_seq .= substr($locus_seq, $block_start[$i], $block_size[$i]);}
			
			# 相補鎖に変換
			&common::complementary($query_seq) if $col[5] < 0;
			
			# クエリー配列を相同性検索に入力
			print SEARCH_IN ">$col[3]\n$query_seq\n";
		}
		
		# ゲノム配列のfastaファイルを閉じる
		close(GENOME);
		
		# 相同性検索への入力を閉じる (-d指定時)
		close(SEARCH_IN) if $opt{"d"};
		
		# パイプを閉じる
		close($report);
		
		# プロセスを終了
		exit(0);
	}
	## ここまで子プロセスの処理 ##
	
	# 相同性検索への入力を閉じる (-d指定時)
	close(SEARCH_IN) if $opt{"d"};
	
	# 変数を宣言
	my %loci = ();
	my %blast_hits = ();
	my @bed_data = (undef, undef, undef, "");
	my $last_query = undef;
	my $last_contig = "";
		
	# パイプを開く
	$report->reader;
	
	# 相同性検索を利用する場合 (-bまたは-d指定時)
	if ($opt{"b"} or $opt{"d"}) {
		# 相同性検索の出力を読み込みながら処理
		print STDERR $opt{"b"} ? "Loading homology search results..." : "Running homology search...";
		while (my $line = <SEARCH_OUT>) {
			# コメント行を除外
			next if $line =~ /^["#"]/;
			
			# 改行コードを除去
			chomp($line);
			
			# タブ文字でデータを分割
			my @col = split(/\t/, $line);
			
			# クエリー名を編集
			($col[0]) = split(/::/, $col[0]);
			
			# 詳細なサブジェクトを追加
			$col[12] = $col[1] if !defined($col[12]);
			
			# サブジェクト名を編集
			($col[1]) = split(/\s/, $col[1]);
			
			# クエリー名が変わった場合
			if (!defined($last_query) or $col[0] ne $last_query) {
				# 条件を満たす場合はハッシュにデータをバイナリ形式で登録
				push(@{$loci{$last_contig}}, pack($bed12_template, @bed_data)) if @{&common::sift_hits(\%blast_hits, $opt{"r"}, $opt{"k"})};
				
				# ヒットハッシュをリセット
				%blast_hits = ();
				
				# クエリー名を更新
				$last_query = $col[0];
				
				# 子プロセスからbedデータを受信しながら処理
				while (my $line = <$report>) {
					# 改行コードを除去
					chomp($line);
					
					# タブ文字でデータを分割
					@bed_data = split(/\t/, $line);
					
					# クエリー名が一致した場合
					last if $bed_data[3] eq $last_query;
				}
				
				# クエリー名が一致しなかった場合
				&exception::error("query not found in bed data probably due to inconsistent data order") if $bed_data[3] ne $last_query;
				
				# データを変換
				$bed_data[10] = pack("L*", split(/,/, $bed_data[10]));
				$bed_data[11] = pack("L*", split(/,/, $bed_data[11]));
				
				# コンティグ名が変わった場合はコンティグ名を更新
				$last_contig = $bed_data[0] if $bed_data[0] ne $last_contig;
			}
			
			# 逆鎖のヒットを除去
			next if ($col[6] < $col[7]) ^ ($col[8] < $col[9]);
			
			# 検索結果をヒットハッシュに登録
			push(@{$blast_hits{$col[1]}}, {"query_start" => $col[6], "query_end" => $col[7], "subject_start" => $col[8], "subject_end" => $col[9], "score" => $col[11], "title" => $col[12]});
		}
		
		# 条件を満たす場合はハッシュに残りのデータをバイナリ形式で登録
		push(@{$loci{$last_contig}}, pack($bed12_template, @bed_data)) if @{&common::sift_hits(\%blast_hits, $opt{"r"}, $opt{"k"})};
		
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
	}
	
	# 相同性検索を利用しない場合 (-bまたは-d未指定時)
	else {
		# 子プロセスからデータを受信しながら処理
		print STDERR "Loading bed data...";
		while (my $line = <$report>) {
			# 改行コードを除去
			chomp($line);
			
			# タブ文字でデータを分割
			my @bed_data = split(/\t/, $line);
			
			# データを変換
			$bed_data[10] = pack("L*", split(/,/, $bed_data[10]));
			$bed_data[11] = pack("L*", split(/,/, $bed_data[11]));
			
			# コンティグ名が変わった場合はコンティグ名を更新
			$last_contig = $bed_data[0] if $bed_data[0] ne $last_contig;
			
			# ハッシュにデータをバイナリ形式で登録
			push(@{$loci{$last_contig}}, pack($bed12_template, @bed_data));
		}
	}
	
	# パイプを閉じる
	close($report);
	
	# 子プロセスが終了するまで待機
	waitpid($pid, 0);
	
	# 子プロセスが異常終了した場合
	&exception::error("process abnormally exited") if $?;
	print STDERR "completed\n";
	
	# 相同性検索でヒットが得られなかった場合 (-d指定時)
	&exception::error("no hits found from $opt{h} search") if $opt{"d"} and !defined($last_query);
	
	# 入出力のキューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入力キューにデータをバイナリ形式で追加
	print STDERR "Running isoform definition...";
	foreach my $contig (keys(%loci)) {$input->enqueue(pack("S/A*L(L/a*)*", $contig, scalar(@{$loci{$contig}}), @{$loci{$contig}})) if @{$loci{$contig}};}
	print STDERR "\rRunning isoform definition...0%";
	
	# 変数を宣言
	my $num_data = $input->pending;
	
	# 指定したワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		threads::async {
			# ワーカースレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューからデータをバイナリ形式で取得して処理
			while (defined(my $dat = $input->dequeue_nb)) {
				# データを変換
				my ($subject, @locus) = unpack("S/A*L/(L/a*)*", $dat);
				
				# アイソフォームを分離
				my $isoforms = &common::define_isoforms(\@locus, $opt{"v"}, $opt{"t"});
				
				# アイソフォームを分離し、出力キューにデータをバイナリ形式で追加
				$output->enqueue(pack("S/A*L(L/a*)*", $subject, scalar(@{$isoforms}), @{$isoforms}));
				
				# 完了データ数を取得
				my $fin_data = $output->pending;
				
				# 経過を表示
				print STDERR "\rRunning isoform definition...", int($fin_data / $num_data * 100), "%" if int($fin_data / $num_data * 100) > int(($fin_data - 1) / $num_data * 100);
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 変数を宣言
	my $fin_threads = 0;
	
	# 各ワーカースレッドが終了するまで待機
	foreach (threads->list) {$_->join and $fin_threads++;}
	
	# ワーカースレッドが異常終了した場合
	&exception::error("thread abnormally exited") if $fin_threads < $opt{"p"};
	
	# 変数を宣言
	my $gene_id = 0;
	
	# 出力キューからデータをバイナリ形式で取得して処理
	while (my $dat = $output->dequeue_nb) {
		# データを変換
		my ($subject, @locus) = unpack("S/A*L/(L/a*)*", $dat);
		
		# 領域データリストにデータをバイナリ形式で追加
		$loci{$subject} = \@locus;
	}
	
	# 結果を出力
	foreach my $subject (sort {&common::faidx_decode($genome_faidx->{$a})->{"id_order"} <=> &common::faidx_decode($genome_faidx->{$b})->{"id_order"}} keys(%loci)) {&common::output($loci{$subject}, $gene_id, $opt{"n"}, $opt{"f"});}
	print STDERR "\rRunning isoform definition...completed\n";
	return(1);
}

## ここからpredictコマンドのパッケージ ##
package predict;

# コマンドとオプションを定義
sub define {
	$note = "Predict protein-coding regions from cDNA sequences under specified conditions.";
	$usage = "<protein.faa> <STDIN|in1.fna> [in2.fna ...] [>out.bed|>out.gtf]";
	$option{"s"} = "\tOmit reverse strand search";
	$option{"g STR "} = "Gene prediction engine <exonerate|genewise> [exonerate]";
	$option{"h STR "} = "Homology search engine <blastx|blastx-fast> [blastx]";
	$option{"k STR "} = "Keywords for filtering (AND[&], OR[;], BUT[!])";
	$option{"l INT "} = "Minimum length to regard as complete CDS <0-> [0]";
	$option{"r INT "} = "Cutoff rank of hits <1->";
	$option{"c FLOAT "} = "Minimum query coverage to regard as complete CDS <0-1> [0.85]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	&exception::error("unknown engine specified: -g $opt{g}") if $opt{"g"} and $opt{"g"} ne "exonerate" and $opt{"g"} ne "genewise";
	&exception::error("unknown engine specified: -h $opt{h}") if $opt{"h"} ne "blastx" and $opt{"h"} ne "blastx-fast";
	&exception::error("specify INT >= 0: -l $opt{l}") if $opt{"l"} !~ /^\d+$/;
	&exception::error("specify INT >= 1: -r $opt{r}") if defined($opt{"r"}) and ($opt{"r"} !~ /^\d+$/ or $opt{"r"} < 1);
	&exception::error("specify FLOAT 0-1: -c $opt{c}") if $opt{"c"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"c"} > 1;
	
	# 依存関係を確認
	my $software_list = &common::check_dependencies(["makeblastdb", "blastx", $opt{"g"}]);
	&exception::error("makeblastdb not installed") if !$software_list->{"makeblastdb"};
	&exception::error("blastx not installed") if !$software_list->{"blastx"};
	&exception::error("$opt{g} not installed") if $opt{"g"} and !$software_list->{$opt{"g"}};
	
	# オプションの対応を確認
	&exception::error("$opt{h} unavailable on your blast+: -h $opt{h}") if $opt{"h"} eq "blastx-fast" and $software_list->{"blastx"} =~ /blast 2\.2\.[1-2]\d/;
	
	# 参照配列のfastaファイル名を取得
	my $ref_file = shift(@ARGV);
	
	# 入力ファイルを確認
	&exception::error("input file not specified") if !@ARGV and !-p STDIN;
	foreach (@ARGV) {
		&exception::error("file not found: $_") if !-f $_;
		&exception::error("file unreadable: $_") if !-r $_;
		&exception::error("null file specified: $_") if !-s $_;
	}
	
	# 参照配列のfastaファイルを検索
	&common::find_db($ref_file);
	
	# 参照配列のfastaインデックスを取得
	my $ref_faidx = &common::read_fasta($ref_file, $opt{"w"});
	
	# 参照配列のblastデータベースを確認 (-b未指定時)
	&common::check_blastdb($ref_file, 1) if !$opt{"b"};
	
	# 検索方向を定義
	my $strand = $opt{"s"} ? "plus" : "both";
	
	# デフォルトの処理を定義
	my %homology_search = (
		"blastx"       => "blastx -num_threads $opt{p} -db $ref_file -outfmt '6 std salltitles' -soft_masking true -strand $strand",
		"blastx-fast"  => "blastx -num_threads $opt{p} -db $ref_file -outfmt '6 std salltitles' -soft_masking true -strand $strand -task blastx-fast"
	);
	my %gene_prediction = (
		"exonerate" => "exonerate -m protein2genome -V 0 --showalignment F --showvulgar F --showtargetgff T",
		"genewise"  => "genewise -divide '#' -pseudo -sum -gff -silent -quiet -nosplice_gtag"
	);
	
	# オプションの処理を追加
	foreach (keys(%homology_search)) {$homology_search{$_} .= `which tee` ? " 2>/dev/null | tee fate_predict_$opt{h}.out" : "";}
	$gene_prediction{"exonerate"} .= $opt{"s"} ? " -r F" : " -r T";
	$gene_prediction{"genewise"} .= $opt{"s"} ? " -tfor" : " -both";
	
	# 個々の配列データを保存しておくディレクトリを作成
	mkdir("prot") or &exception::error("failed to make directory: prot") if !-d "prot";
	mkdir("nucl") or &exception::error("failed to make directory: nucl") if !-d "nucl";
	
	# 相同性検索の出力ファイルを開く (-b指定時)
	open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}") if $opt{"b"};
	
	# 相同性検索を実行 (-b未指定時)
	main::open2(*SEARCH_OUT, *SEARCH_IN, "$homology_search{$opt{h}} 2>/dev/null") or &exception::error("failed to execute homology search") if !$opt{"b"};
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	&exception::error("failed to fork process") if !defined($pid);
	
	## ここから子プロセスの処理 ##
	if (!$pid) {
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# 変数を宣言
		my $query_title = undef;
		my $query_seq = "";
		
		# 入力の改行コードを一時的に変更 (-w指定時)
		local $/ = "\r\n" if $opt{"w"};
		
		# クエリー配列を読み込みながら処理
		while (my $line = <>) {
			# 改行コードを除去
			chomp($line);
			
			# 配列行の処理
			$query_seq .= $line =~ /^>/ ? "" : uc($line);
			
			# ID行の処理
			if ($line =~ /^>/ or eof) {
				# 配列データを処理
				if (defined($query_title)) {
					# 出力ファイルを作成
					open(NUCL, ">", "nucl/$query_title.fa") or &exception::error("failed to make file: nucl/$query_title.fa");
					
					# クエリー配列をfasta形式でファイルに出力
					print NUCL ">$query_title\n$query_seq\n";
					
					# 出力ファイルを閉じる
					close(NUCL);
					
					# クエリー配列を相同性検索に入力 (-b未指定時)
					print SEARCH_IN ">$query_title\n$query_seq\n" if !$opt{"b"};
				}
				
				# ID行の最初の空白文字の前までを配列名として登録
				($query_title) = split(/\s/, substr($line, 1));
				
				# 配列をリセット
				$query_seq = "";
			}
		}
		
		# 相同性検索への入力を閉じる (-b未指定時)
		close(SEARCH_IN) if !$opt{"b"};
		
		# プロセスを終了
		exit(0);
	}
	## ここまで子プロセスの処理 ##
	
	# 相同性検索への入力を閉じる (-b未指定時)
	close(SEARCH_IN) if !$opt{"b"};
	
	# 変数を宣言
	my %filtered_hits = ();
	my %blast_hits = ();
	my $last_query = undef;
	
	# 相同性検索の出力を読み込みながら処理
	print STDERR $opt{"b"} ? "Loading homology search results..." : "Running homology search...";
	while (my $line = <SEARCH_OUT>) {
		# コメント行を除外
		next if $line =~ /^["#"]/;
		
		# 改行コードを除去
		chomp($line);
		
		# タブ文字でデータを分割
		my @col = split(/\t/, $line);
		
		# 詳細なサブジェクトを追加
		$col[12] = $col[1] if !defined($col[12]);
		
		# サブジェクト名を編集
		($col[1]) = split(/\s/, $col[1]);
		
		# クエリー名が変わった場合
		if (!defined($last_query) or $col[0] ne $last_query) {
			# 条件を満たす各ヒットについてハッシュにデータを登録
			foreach my $subject (@{&common::sift_hits(\%blast_hits, $opt{"r"}, $opt{"k"})}) {$filtered_hits{$last_query}{$subject} = $blast_hits{$subject};}
			
			# ヒットハッシュをリセット
			%blast_hits = ();
			
			# クエリー名を更新
			$last_query = $col[0];
		}
		
		# 検索結果をヒットハッシュに登録
		push(@{$blast_hits{$col[1]}}, {"query_start" => $col[6], "query_end" => $col[7], "subject_start" => $col[8], "subject_end" => $col[9], "score" => $col[11], "title" => $col[12]});
	}
	
	# 条件を満たす各ヒットについてハッシュにデータを登録
	foreach my $subject (@{&common::sift_hits(\%blast_hits, $opt{"r"}, $opt{"k"})}) {$filtered_hits{$last_query}{$subject} = $blast_hits{$subject};}
	
	# 相同性検索の出力を閉じる
	close(SEARCH_OUT);
	
	# 子プロセスが終了するまで待機
	waitpid($pid, 0);
	
	# 子プロセスが異常終了した場合
	&exception::error("process abnormally exited") if $?;
	print STDERR "completed\n";
	
	# 相同性検索でヒットが得られなかった場合
	&exception::error("no hits found from $opt{h} search") if !defined($last_query);
	
	# 変数を宣言
	my %loci = ();
	
	# 入出力のキューを作成
	my $input = Thread::Queue->new;
	my $output = Thread::Queue->new;
	
	# 入力キューにデータをバイナリ形式で追加
	print STDERR "Running gene prediction...";
	foreach my $query (keys(%filtered_hits)) {
		foreach my $subject (keys(%{$filtered_hits{$query}})) {
			$input->enqueue(pack("S/A*S/A*C", $query, $subject, List::Util::reduce {$a | $b} map {(($_->{"query_start"} < $_->{"query_end"}) ^ ($_->{"subject_start"} < $_->{"subject_end"})) + 1} @{$filtered_hits{$query}{$subject}}));
		}
	}
	print STDERR "\rRunning gene prediction...0%";
	
	# 変数を宣言
	my $num_data = $input->pending;
	
	# 指定したワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		threads::async {
			# ワーカースレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 参照配列のfastaファイルを開く
			open(REF, "<", $ref_file) or &exception::error("failed to open file: $ref_file");
			
			# 入力キューからデータをバイナリ形式で取得して処理
			while (defined(my $dat = $input->dequeue_nb)) {
				# データを変換
				my ($query, $subject, $strand) = unpack("S/A*S/A*C", $dat);
				
				# 変数を宣言
				my $target_seq = "";
				
				# 対象配列ファイルを開く
				open(TARGET, "<", "nucl/$query.fa") or &exception::error("failed to open file: nucl/$query.fa");
				
				# 対象配列ファイルを読み込みながら処理
				while (my $line = <TARGET>) {
					# 改行コードを除く
					chomp($line);
					
					# 配列を取得
					$target_seq .= $line if $line !~ /^>/;
				}
				
				# 対象配列ファイルを閉じる
				close(TARGET);
				
				# 変数を宣言
				my $faidx = &common::faidx_decode($ref_faidx->{$subject});
				
				# 参照配列を取得
				seek(REF, $faidx->{"seq_start"}, 0);
				read(REF, my $ref_seq, int($faidx->{"seq_length"} / $faidx->{"row_width"}) * $faidx->{"row_bytes"} + $faidx->{"seq_length"} % $faidx->{"row_width"});
				
				# 配列中に含まれる改行コードとアスタリスクを除去
				$ref_seq =~ s/\n|\r|\*//g;
				
				# 配列を大文字に変換
				$ref_seq = uc($ref_seq);
				
				# 参照配列ファイルを作成
				open(PROT, ">", "prot/ref$thread_id.fa") or &exception::error("failed to make file: prot/ref$thread_id.fa");
				
				# 取得した配列をファイルに出力
				print PROT ">$subject\n$ref_seq\n";
				
				# 参照配列ファイルを閉じる
				close(PROT);
				
				# 変数を宣言
				my @genes = ();
				my $query_file = "";
				my $target_file = "";
				my $query_start = 0;
				my $query_end = 0;
				my $summary_flag = 0;
				
				# 遺伝子構造予測の引数を定義
				$query_file = "-q prot/ref$thread_id.fa" if $opt{"g"} eq "exonerate";
				$query_file = "prot/ref$thread_id.fa" if $opt{"g"} eq "genewise";
				$target_file = "-t nucl/$query.fa" if $opt{"g"} eq "exonerate";
				$target_file = "nucl/$query.fa" if $opt{"g"} eq "genewise";
				$target_file =~ s/\|/\\\|/g;
				
				# 遺伝子構造予測を実行
				open(PREDICT, "-|", "$gene_prediction{$opt{g}} $query_file $target_file 2>/dev/null") or &exception::error("failed to execute gene prediction: $subject vs $query");
				
				# 遺伝子構造予測結果を1行ずつ読み込んで処理
				while (my $line = <PREDICT>) {
					# コメント行を除外
					++$summary_flag and next if $line =~ /^["#"]/;
					
					# 空白文字でデータを分割
					my @col = split(/\s+/, $line, 9);
					
					# 見出し行を除外
					next if $col[0] eq "Bits" and $col[1] eq "Query" and $col[4] eq "Target";
					
					# summary行の処理 (-g genewise指定時)
					($query_start, $query_end) = (($col[2] - 1) * 3, $col[3] * 3) and next if $opt{"g"} eq "genewise" and !($summary_flag & 0x01);
					
					# 相同性検索がマッチしなかった向きの予測を除外
					next if !($strand & ($col[6] eq "-") + 1);
					
					# 開始点 > 終止点の場合は開始点と終止点を入れ替える
					($col[3], $col[4]) = ($col[4], $col[3]) if $col[3] > $col[4];
					
					# 開始点をbed形式に合わせる
					$col[3]--;
					
					# match行またはgene行の処理
					push(@genes, [$query, $col[3], $col[4], $subject, $col[5], $col[6], 0, 0, ".", 0, [], [], $query_start, $query_end]) if $col[2] eq "match" or $col[2] eq "gene";
					
					# cds行の処理
					++$genes[-1]->[9] and push(@{$genes[-1]->[10]}, $col[4] - $col[3]) and push(@{$genes[-1]->[11]}, $col[3]) if $col[2] eq "cds";
					
					# similarity行の処理 (-g exonerate指定時)
					my @query_block = map {[split(/\s/)]} grep {/^Align/} split(/ ; /, $col[8]) if $opt{"g"} eq "exonerate" and $col[2] eq "similarity";
					$genes[-1]->[12] = List::Util::min(map {($_->[2] - 1) * 3} @query_block) if @query_block;
					$genes[-1]->[13] = List::Util::max(map {($_->[2] - 1) * 3 + $_->[3]} @query_block) if @query_block;
				}
				
				# 遺伝子構造予測を終了
				close(PREDICT);
				
				# 各遺伝子について翻訳領域を推定
				foreach my $gene (@genes) {
					# 変数を宣言
					my $sign = $gene->[5] eq "-";
					my $completeness = 0;
					
					# 上流配列を取得
					my $upstream_seq = substr($target_seq, 0, $gene->[1]);
					
					# 下流配列を取得
					my $downstream_seq = substr($target_seq, $gene->[2]);
					
					# 逆向きの場合は上流配列と下流配列を入れ替え
					($upstream_seq, $downstream_seq) = ($downstream_seq, $upstream_seq) if $sign;
					
					# 5'側クエリー被覆率が指定値以上の場合
					if ($faidx->{"seq_length"} - $gene->[12] >= $faidx->{"seq_length"} * $opt{"c"}) {
						# 上流フランキング領域のアミノ酸配列を取得
						my $upstream_aa = &common::translate($upstream_seq, $gene->[5] . length($upstream_seq) % 3 - $sign);
						
						# 先頭エキソンのアミノ酸配列を取得
						my $first_exon_aa = &common::translate(substr($target_seq, $gene->[11]->[0], $gene->[10]->[0]), -$sign);
						
						# 開始コドンを探索
						my $start_pos = rindex($upstream_aa, "M");
						my $stop_pos = rindex($upstream_aa, "*");
						$start_pos = $stop_pos < $start_pos ? $start_pos - length($upstream_aa) : index($first_exon_aa, "M") + 1;
						$start_pos -= $start_pos > 0;
						$start_pos *= 3;
						
						# 開始コドンを考慮しても5'側クエリー被覆率が指定値以上の場合は先頭ブロックを修正
						($gene->[10]->[0], $gene->[11]->[0]) = ($gene->[10]->[0] - $start_pos * (1 - $sign), $gene->[11]->[0] + $start_pos) if $faidx->{"seq_length"} - $gene->[12] - $start_pos >= $faidx->{"seq_length"} * $opt{"c"};
						
						# 5'端が完全であることを認定
						$completeness++;
					}
					
					# 3'側クエリー被覆率が指定値以上の場合
					if ($gene->[13] >= $faidx->{"seq_length"} * $opt{"c"}) {
						# 下流フランキング領域のアミノ酸配列を取得
						my $downstream_aa = &common::translate($downstream_seq, -$sign);
						
						# 末尾エキソンのアミノ酸配列を取得
						my $last_exon_aa = &common::translate(substr($target_seq, $gene->[11]->[-1], $gene->[10]->[-1]), $gene->[5] . $gene->[10]->[-1] % 3 - $sign);
						
						# 終止コドンを探索
						my $terminal_pos = rindex($last_exon_aa, "*");
						$terminal_pos = $terminal_pos < 0 ? index($downstream_aa, "*") : $terminal_pos - length($last_exon_aa);
						$terminal_pos++;
						$terminal_pos *= 3;
						
						# 終止コドンを考慮しても3'側クエリー被覆率が指定値以上の場合は末尾ブロックを修正
						($gene->[10]->[-1], $gene->[11]->[-1]) = ($gene->[10]->[-1] + $terminal_pos, $gene->[11]->[-1] - $terminal_pos * $sign) if $gene->[13] + $terminal_pos >= $faidx->{"seq_length"} * $opt{"c"};
						
						# 3'端が完全であることを認定
						$completeness++;
					}
					
					# 逆向きの場合はブロックを並べ替え
					($gene->[10], $gene->[11]) = ([reverse(@{$gene->[10]})], [reverse(@{$gene->[11]})]) if $sign;
					
					# 領域開始点・終了点を修正
					($gene->[1], $gene->[2]) = ($gene->[11]->[0], $gene->[10]->[-1] + $gene->[11]->[-1]);
					
					# ブロックの基準点を修正
					my $basal_pos = $gene->[11]->[0];
					
					# ブロックの相対位置を修正
					$gene->[11] = [map {$_ - $basal_pos} @{$gene->[11]}];
					
					# 翻訳領域のアミノ酸配列を取得
					my $cds = "";
					for (my $i = 0;$i < $gene->[9];$i++) {$cds .= substr($target_seq, $gene->[1] + $gene->[11]->[$i], $gene->[10]->[$i]);}
					my $aa = &common::translate($cds, -$sign);
					
					# フレームシフトが存在する場合 (偽遺伝子)
					if (List::Util::sum(@{$gene->[10]}) % 3 > 0) {$gene->[8] = "red";}
					
					# 終止コドンが末尾以外に存在する場合 (偽遺伝子)
					elsif (index($aa, "*") >= 0 and index($aa, "*") < length($aa) - 1) {$gene->[8] = "red";}
					
					# 両末端が完全とみなされない場合 (偽遺伝子)
					elsif ($completeness < 2) {$gene->[8] = "red";}
					
					# 配列長が指定値未満の場合 (偽遺伝子)
					elsif (List::Util::sum(@{$gene->[10]}) < $opt{"l"}) {$gene->[8] = "red";}
					
					# 上記に該当しない場合 (機能遺伝子)
					else {$gene->[8] = "blue";}
					
					# データを修正
					$gene->[5] .= 1;
					($gene->[6], $gene->[7]) = ($gene->[1], $gene->[2]);
					
					# データをバイナリ形式に変換
					$gene = pack($bed12_template, @{$gene}[0..9], pack("L*", @{$gene->[10]}), pack("L*", @{$gene->[11]}));
				}
				
				# 出力キューにデータをバイナリ形式で追加
				$output->enqueue(pack("S/A*L(L/a*)*", $query, scalar(@genes), @genes));
				
				# 完了データ数を取得
				my $fin_data = $output->pending;
				
				# 経過を表示
				print STDERR "\rRunning gene prediction...", int($fin_data / $num_data * 100), "%" if int($fin_data / $num_data * 100) > int(($fin_data - 1) / $num_data * 100);
			}
			
			# 参照配列のfastaファイルを閉じる
			close(REF);
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 変数を宣言
	my $fin_threads = 0;
	
	# 各ワーカースレッドが終了するまで待機
	foreach (threads->list) {$_->join and $fin_threads++;}
	
	# ワーカースレッドが異常終了した場合
	&exception::error("thread abnormally exited") if $fin_threads < $opt{"p"};
	
	# 出力キューからデータをバイナリ形式で取得して処理
	while (my $dat = $output->dequeue_nb) {
		# データを変換
		my ($query, @locus) = unpack("S/A*L/(L/a*)*", $dat);
		
		# ハッシュにデータをバイナリ形式で登録
		push(@{$loci{$query}}, @locus) if @locus;
	}
	print STDERR "\rRunning gene prediction...completed\n";
	
	# 入力キューにアセンブルデータをバイナリ形式で追加
	print STDERR "Running isoform definition...";
	foreach my $query (keys(%loci)) {$input->enqueue(pack("S/A*L(L/a*)*", $query, scalar(@{$loci{$query}}), @{$loci{$query}}));}
	print STDERR "\rRunning isoform definition...0%";
	
	# 変数をリセット
	$num_data = $input->pending;
	
	# 指定したワーカースレッド数で並列処理
	for (my $thread_id = 0;$thread_id < $opt{"p"};$thread_id++) {
		## ここからワーカースレッドの処理 ##
		threads::async {
			# ワーカースレッドのみを終了可能に変更
			threads->set_thread_exit_only(1);
			
			# 入力キューからデータをバイナリ形式で取得して処理
			while (defined(my $dat = $input->dequeue_nb)) {
				# データを変換
				my ($query, @locus) = unpack("S/A*L/(L/a*)*", $dat);
				
				# アイソフォームを分離
				my $isoforms = &common::define_isoforms(\@locus, $opt{"v"}, $opt{"t"});
				
				# アイソフォームを分離し、出力キューにデータをバイナリ形式で追加
				$output->enqueue(pack("S/A*L(L/a*)*", $query, scalar(@{$isoforms}), @{$isoforms}));
				
				# 完了データ数を取得
				my $fin_data = $output->pending;
				
				# 経過を表示
				print STDERR "\rRunning isoform definition...", int($fin_data / $num_data * 100), "%" if int($fin_data / $num_data * 100) > int(($fin_data - 1) / $num_data * 100);
			}
			
			# スレッドを終了
			return(1);
		};
		## ここまでワーカースレッドの処理 ##
	}
	
	# 変数をリセット
	$fin_threads = 0;
	
	# 各ワーカースレッドが終了するまで待機
	foreach (threads->list) {$_->join and $fin_threads++;}
	
	# ワーカースレッドが異常終了した場合
	&exception::error("thread abnormally exited") if $fin_threads < $opt{"p"};
	
	# 変数を宣言
	my $gene_id = 0;
	
	# 出力キューからデータをバイナリ形式で取得して処理
	while (my $dat = $output->dequeue_nb) {
		# データを変換
		my ($query, @locus) = unpack("S/A*L/(L/a*)*", $dat);
		
		# 領域データリストにデータをバイナリ形式で追加
		$loci{$query} = \@locus;
	}
	
	# 結果を出力
	foreach my $query (sort {$a cmp $b} keys(%loci)) {&common::output($loci{$query}, $gene_id, $opt{"n"}, $opt{"f"});}
	print STDERR "\rRunning isoform definition...completed\n";
	return(1);
}

## ここから共通処理のパッケージ ##
package common;

# キーワード検索 common::keyword_search(検索対象文字列, キーワード)
sub keyword_search {
	# キーワード未指定時はNULLを返す
	return("") if !$_[1];
	
	# 検索を実行して結果を返す
	return(List::Util::first {!defined(List::Util::first {if (/^!/) {s/^!//;$_[0] =~ /$_/i} else {$_[0] !~ /$_/i}} split(/&/, $_))} split(/;/, $_[1]));
}

# 相補鎖変換 common::complementary(配列)
sub complementary {
	# 配列を逆順に並べ替える
	$_[0] = reverse($_[0]);
	
	# 相補的な塩基に置換
	$_[0] =~ tr/ATGCRYKMDBVH/TACGYRMKHVBD/;
	return(1);
}

# アミノ酸変換 common::translate(配列, 読み枠)
sub translate {
	# 変数を宣言
	my $aa = "";
	
	# 読み枠が負の値の場合は配列を相補鎖変換
	&common::complementary($_[0]) if $_[1] < 0;
	
	# 3塩基ごとにアミノ酸配列に変換
	for (my $site = abs($_[1]) - ($_[1] < 0);$site < length($_[0]);$site += 3) {
		my $triplet = substr($_[0], $site, 3);
		$aa .= exists($codon{$triplet}) ? $codon{$triplet} : "X";
	}
	
	# 翻訳したアミノ酸配列を返す
	return($aa);
}

# 立っているビット数を算出 common::population_count(ビット列)
sub population_count {
	# 変数を宣言
	my $count = 0;
	
	# 64-bitごとに分割統治法で立っているビット数を算出
	foreach (unpack("Q*", $_[0] . chr(0) x 7)) {
		$_ = ($_ & 0x5555555555555555) + ($_ >> 1 & 0x5555555555555555);
		$_ = ($_ & 0x3333333333333333) + ($_ >> 2 & 0x3333333333333333);
		$_ = ($_ & 0x0F0F0F0F0F0F0F0F) + ($_ >> 4 & 0x0F0F0F0F0F0F0F0F);
		$_ = ($_ & 0x00FF00FF00FF00FF) + ($_ >> 8 & 0x00FF00FF00FF00FF);
		$_ = ($_ & 0x0000FFFF0000FFFF) + ($_ >> 16 & 0x0000FFFF0000FFFF);
		$_ = ($_ & 0x00000000FFFFFFFF) + ($_ >> 32 & 0x00000000FFFFFFFF);
		$count += $_;
	}
	
	# 結果を返す
	return($count);
}

# データベース検索 common::find_db(ファイル名)
sub find_db {
	# 環境変数からデータベースパスを取得
	my @db_list = $ENV{"BLASTDB"} ? map {"$_/"} split(/:/, $ENV{"BLASTDB"}) : ();
	
	# ファイルを探索
	my $prefix = List::Util::first {-f "$_$_[0]"} ("", @db_list);
	
	# ファイルが存在しない場合
	&exception::error("file not found: $_[0]") if !defined($prefix);
	
	# ファイルパスを変更
	$_[0] = "$prefix$_[0]";
	
	# ファイルを確認
	&exception::error("file unreadable: $_[0]") if !-r $_[0];
	&exception::error("null file specified: $_[0]") if !-s $_[0];
	return(1);
}

# BLASTデータベース確認 common::check_blastdb(ファイルパス, 配列型)
sub check_blastdb {
	# 配列型と検索コマンドを定義
	my ($seq_type, $command) = $_[1] ? ("prot", "blastx") : ("nucl", "blastn");
	
	# BLASTデータベースファイルを確認
	return(1) if !system("$command -query /dev/null -db $_[0] >/dev/null 2>&1");
	&exception::caution("BLAST database file unavailable: $_[0]");
	
	# BLASTデータベースファイルを作成
	print STDERR "Building BLAST database...";
	&exception::error("failed to build BLAST database") if system("makeblastdb -in $_[0] -dbtype $seq_type 1>/dev/null 2>&1");
	print STDERR "completed\n";
	return(1);
}

# 依存関係確認 common::check_dependencies(コマンドリストリファレンス)
sub check_dependencies {
	# 変数を宣言
	my %result = ();
	
	# 各コマンドについて処理
	foreach (@{$_[0]}) {
		next if !$_;
		$result{$_} .= `$_ 2>/dev/null`;
		$result{$_} .= `$_ -version 2>/dev/null`;
	}
	
	# 結果を返す
	return(\%result);
}

# bed形式確認 common::check_bed(bedデータリストリファレンス, 領域名プレフィックス)
sub check_bed {
	# 3列未満のデータの場合
	&exception::error("bad format input") if @{$_[0]} < 3;
	
	# 4列目が未定義の場合
	$_[0]->[3] = $_[1] if !defined($_[0]->[3]);
	
	# 5列目が未定義の場合
	$_[0]->[4] = 1000 if !defined($_[0]->[4]);
	
	# 6列目が未定義の場合
	$_[0]->[5] = "+" if !defined($_[0]->[5]);
	
	# 7列目が未定義の場合
	$_[0]->[6] = $_[0]->[1] if !defined($_[0]->[6]);
	
	# 8列目が未定義の場合
	$_[0]->[7] = $_[0]->[2] if !defined($_[0]->[7]);
	
	# 9列目が未定義の場合
	$_[0]->[8] = "." if !defined($_[0]->[8]);
	
	# 10列目が未定義の場合
	$_[0]->[9] = 1 if !defined($_[0]->[9]);
	
	# 11列目が未定義の場合
	$_[0]->[10] = $_[0]->[2] - $_[0]->[1] if !defined($_[0]->[10]);
	
	# 12列目が未定義の場合
	$_[0]->[11] = 0 if !defined($_[0]->[11]);
	
	# 10列目の値と11列目および12列目の要素数が一致していない場合
	&exception::error("bad format input") if $_[0]->[9] != split(/,/, $_[0]->[10]) or $_[0]->[9] != split(/,/, $_[0]->[11]);
	
	# 向きを数字に変換
	$_[0]->[5] = $_[0]->[5] . "1";
	return(1);
}

# fastaインデックス読み込み common::read_fasta(ファイルパス, 改行コード型)
sub read_fasta {
	# 変数を宣言
	my %fasta_index = ();
	my $id_key = "";
	my $id_order = 0;
	my $id_start = 0;
	my $seq_length = 0;
	my $seq_start = 0;
	my $row_width = 0;
	my $row_bytes = 0;
	my $file_check = 0;
	
	# fastaインデックスファイルを確認
	$file_check = &exception::caution("index file not found: $_[0].fai") if !$file_check and !-f "$_[0].fai";
	$file_check = &exception::caution("index file unreadable: $_[0].fai") if !$file_check and !-r "$_[0].fai";
	$file_check = &exception::caution("null index file specified: $_[0].fai") if !$file_check and !-s "$_[0].fai";
	
	# インデックスファイルが存在する場合
	if (!$file_check) {
		# fastaインデックスファイルを開く
		open(FASTA_INDEX, "<", "$_[0].fai") or &exception::error("failed to open index file: $_[0].fai");
		
		# データを読み込みながら処理
		print STDERR "Loading fasta index...";
		while (my $line = <FASTA_INDEX>) {
			# 改行コードを除去
			chomp($line);
			
			# タブ文字でデータを分割
			($id_key, $seq_length, $seq_start, $row_width, $row_bytes) = split(/\t/, $line);
			
			# インデックスハッシュにデータをバイナリ形式で登録
			$fasta_index{$id_key} = pack("QQLQLL", $id_order, $id_start, $seq_length, $seq_start, $row_width, $row_bytes);
			
			# ID開始点を更新
			$id_start = $seq_start + int($seq_length / $row_width) * $row_bytes + $seq_length % $row_width;
			
			# 行番号を加算
			$id_order++;
		}
		print STDERR "completed\n";
		
		# fastaインデックスファイルを閉じる
		close(FASTA_INDEX);
	}
	
	# インデックスファイルが存在しない場合
	else {
		# 変数を宣言
		my $error_flag = 0;
		
		# fastaファイルを開く
		open(FASTA, "<", $_[0]) or &exception::error("failed to open file: $_[0]");
		
		# fastaインデックスファイルを新規作成
		open(FASTA_INDEX, ">", "$_[0].fai") or &exception::error("failed to build index file: $_[0].fai");
		
		# 入力の改行コードを一時的に変更 (-w指定時)
		local $/ = "\r\n" if $_[1];
		
		# データを読み込みながら処理
		print STDERR "Building fasta index...";
		while (my $line = <FASTA>) {
			# 1行の文字列数を取得
			my $row_bytes1 = length($line);
			
			# 改行コードを除去
			chomp($line);
			
			# 配列行の処理
			if ($id_key and $line !~ /^>/) {
				# データの整合性を確認
				&exception::error("different line length detected: $id_key", $_[0]) if $error_flag == 1;
				&exception::error("different EOL code detected: $id_key", $_[0]) if $error_flag == 2;
				$error_flag = 1 if $row_width and $row_width != length($line);
				$error_flag = 2 if $row_bytes and $row_bytes != $row_bytes1;
				
				# 配列長を追加
				$seq_length += length($line);
				$row_width = length($line) if !$row_width;
				$row_bytes = $row_bytes1 if !$row_bytes;
			}
			
			# ID行およびファイル末の処理
			if ($line =~ /^>/ or eof(FASTA)) {
				if ($id_key) {
					# インデックスハッシュに直前のデータをバイナリ形式で登録
					$fasta_index{$id_key} = pack("QQLQLL", $id_order, $id_start, $seq_length, $seq_start, $row_width, $row_bytes);
					
					# インデックスファイルに直前のデータを書き込む
					print FASTA_INDEX "$id_key\t$seq_length\t$seq_start\t$row_width\t$row_bytes\n";
				}
				
				# IDの最初の空白文字の直前までをハッシュキーとして登録
				($id_key) = split(/\s/, substr($line, 1));
				
				# 配列開始点を更新
				$seq_start = tell(FASTA);
				
				# ID開始点を更新
				$id_start = $seq_start - $row_bytes1;
				
				# 行番号を加算
				$id_order++;
				
				# その他の情報をリセット
				$seq_length = 0;
				$row_width = 0;
				$row_bytes = 0;
				$error_flag = 0;
			}
		}
		print STDERR "completed\n";
		
		# fastaインデックスファイルを閉じる
		close(FASTA_INDEX);
		
		# fastaファイルを閉じる
		close(FASTA);
	}
	
	# インデックスハッシュのリファレンスを返す
	return(\%fasta_index);
}

# fastaインデックス復号 common::faidx_decode(fastaインデックス)
sub faidx_decode {
	# 変数を宣言
	my %fasta_index = ();
	
	# データを変換してハッシュに登録
	@fasta_index{("id_order", "id_start", "seq_length", "seq_start", "row_width", "row_bytes")} = unpack("QQLQLL", $_[0]);
	
	# ハッシュリファレンスを返す
	return(\%fasta_index);
}

# 相同性検索ヒットの選別 common::sift_hits(相同性検索ヒットハッシュリファレンス, ランク閾値, キーワード)
sub sift_hits {
	# 変数を宣言
	my %score = ();
	
	# 各サブジェクトについて処理
	foreach my $subject (keys(%{$_[0]})) {
		# 変数を宣言
		my $query_array = "";
		my $target_array = "";
		my $total_score = 0;
		my $total_length = 0;
		
		# 各ヒットについて処理
		foreach my $hit (@{$_[0]->{$subject}}) {
			for (List::Util::min($hit->{"query_start"}, $hit->{"query_end"})..List::Util::max($hit->{"query_start"}, $hit->{"query_end"})) {vec($query_array, $_, 1) = 1;}
			for (List::Util::min($hit->{"subject_start"}, $hit->{"subject_end"})..List::Util::max($hit->{"subject_start"}, $hit->{"subject_end"})) {vec($target_array, $_, 1) = 1;}
			$total_score += $hit->{"score"};
			$total_length += abs($hit->{"query_end"} - $hit->{"query_start"}) + abs($hit->{"subject_end"} - $hit->{"subject_start"}) + 2;
		}
		
		# スコアを算出
		$score{$subject} = &common::population_count($query_array . $target_array) * $total_score / $total_length;
	}
	
	# ヒットをスコア順に並べ替え
	my @subjects = sort {$score{$b} <=> $score{$a}} keys(%score);
	
	# 指定値より下位のヒットを削除 (ランク閾値指定時)
	splice(@subjects, $_[1]) if $_[1];
	
	# キーワード条件を満たすサブジェクトリストリファレンスを返す
	return([grep {defined(&common::keyword_search($_[0]->{$_}->[0]->{"title"}, $_[2]))} @subjects]);
}

# アイソフォームを決定 common::define_isoforms(領域データリストリファレンス, 出力アイソフォーム数, 出力遺伝子型)
sub define_isoforms {
	# 変数を宣言
	my %type = ("blue" => $_[2] & 0x01, "yellow" => $_[2] & 0x02, "red" => $_[2] & 0x04, "." => 1);
	my @loci = ();
	
	# データを順に整理
	foreach my $locus (map {[unpack($bed12_template, $_)]} @{$_[0]}) {
		# データを変換
		my @block_size = unpack("L*", $locus->[10]);
		my @block_start = unpack("L*", $locus->[11]);
		
		# 未指定の遺伝子タイプのデータを除外
		next if !$type{$locus->[8]};
		
		# ブロック情報を絶対位置に変換して最後列に追加
		my @blocks = ();
		for (my $i = 0;$i < $locus->[9];$i++) {$blocks[$i] = [$locus->[1] + $block_start[$i], $locus->[1] + $block_size[$i] + $block_start[$i]];}
		
		# 変数を宣言
		my $basal_locus = undef;
		
		# 既出領域と比較
		for (my $i = 0;$i < @loci;$i++) {
			# 方向が異なる場合を除外
			next if $locus->[5] != $loci[$i]->[0]->[3];
			
			# 領域がオーバーラップしない場合を除外
			next if $locus->[2] <= $loci[$i]->[0]->[1] or $loci[$i]->[0]->[2] <= $locus->[1];
			
			# ブロックがオーバーラップしない場合を除外
			next if !grep {my $block = $_;List::Util::first {$_->[0] < $block->[1] and $block->[0] < $_->[1]} @{$loci[$i]->[0]->[4]}} @blocks;
			
			# 既存アイソフォームから完全に一致するものを探索
			my $synonym = List::Util::first {$locus->[1] == $_->[1] and $locus->[10] eq $_->[10] and $locus->[11] eq $_->[11]} @{$loci[$i]->[1]};
			
			# 既存アイソフォームと完全に一致する場合を除外
			if (defined($synonym)) {
				($synonym->[3], $synonym->[4], $synonym->[8]) = ($locus->[3], $locus->[4], $locus->[8]) if $locus->[4] > $synonym->[4] or $locus->[4] == $synonym->[4] and $locus->[3] lt $synonym->[3];
				$basal_locus = 1;
				last;
			}
			
			# オーバーラップする領域を統合
			if (defined($basal_locus)) {
				$basal_locus->[0]->[1] = $loci[$i]->[0]->[1] if $loci[$i]->[0]->[1] < $basal_locus->[0]->[1];
				$basal_locus->[0]->[2] = $loci[$i]->[0]->[2] if $loci[$i]->[0]->[2] > $basal_locus->[0]->[2];
				push(@{$basal_locus->[0]->[4]}, @{$loci[$i]->[0]->[4]});
				push(@{$basal_locus->[1]}, @{$loci[$i]->[1]});
				splice(@loci, $i, 1);
				$i--;
			}
			
			# オーバーラップする領域に追加
			else {
				$basal_locus = $loci[$i];
				$basal_locus->[0]->[1] = $locus->[1] if $locus->[1] < $basal_locus->[0]->[1];
				$basal_locus->[0]->[2] = $locus->[2] if $locus->[2] > $basal_locus->[0]->[2];
				push(@{$basal_locus->[0]->[4]}, @blocks);
				push(@{$basal_locus->[1]}, $locus);
			}
		}
		
		# 新規領域を登録
		push(@loci, [[$locus->[0], $locus->[1], $locus->[2], $locus->[5], \@blocks], [$locus]]) if !defined($basal_locus);
	}
	
	# 各領域について処理
	foreach my $locus (@loci) {
		# 出力アイソフォーム数が指定されている場合
		if ($_[1]) {
			# アイソフォームを並べ替える (スコア > 領域長 > 名前)
			@{$locus->[1]} = sort {$b->[4] <=> $a->[4] or $b->[2] - $b->[1] <=> $a->[2] - $a->[1] or $a->[3] cmp $b->[3]} @{$locus->[1]};
			
			# 上位から指定した個数だけ選抜
			splice(@{$locus->[1]}, $_[1]);
			
			# 領域の位置情報を更新
			$locus->[0]->[1] = List::Util::min(map {$_->[1]} @{$locus->[1]});
			$locus->[0]->[2] = List::Util::max(map {$_->[2]} @{$locus->[1]});
		}
		
		# アイソフォームを並べ替える (開始点 > 終了点 > 名前 > スコア)
		@{$locus->[1]} = sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] or $a->[3] cmp $b->[3] or $b->[4] <=> $a->[4]} @{$locus->[1]};
		
		# 領域の全ブロック情報を削除
		pop(@{$locus->[0]});
	}
	
	# 各領域についてデータをバイナリ形式に変換
	foreach my $locus (@loci) {@{$locus->[1]} = map {pack($bed12_template, @{$_})} @{$locus->[1]};}
	
	# 並べ替え (開始点 > 終了点 > 向き) た領域データリストリファレンスを返す
	return([map {pack("S/A*LLcL(L/a*)*", @{$_->[0]}, scalar(@{$_->[1]}), @{$_->[1]})} sort {$a->[0]->[1] <=> $b->[0]->[1] or $a->[0]->[2] <=> $b->[0]->[2] or $a->[0]->[3] cmp $b->[0]->[3]} @loci]);
}

# データを出力 common::output(領域データリスト, 遺伝子ID, プレフィックス, 出力形式)
sub output {
	# 変数を宣言
	my %type = ("red" => "pseudogene", "yellow" => "truncated", "blue" => "protein_coding", "." => "uncharacterized");
	
	# データを出力
	foreach my $locus (map {[unpack("S/A*LLcL/(L/a*)*", $_)]} @{$_[0]}) {
		# 遺伝子IDを更新
		$_[1]++;
		
		# 遺伝子行を出力 (gtf形式指定時)
		print "$locus->[0]\tfate\tgene\t", $locus->[1] + 1, "\t$locus->[2]\t.\t", $locus->[3] > 0 ? "+" : "-", "\t.\t", 'gene_id "', "$_[2]$_[1]", '";', "\n" if $_[3] eq "gtf";
		
		# 遺伝子行データを削除
		splice(@{$locus}, 0, 4);
		
		# 変数を宣言
		my $isoform_id = 0;
		
		# アイソフォームを順に処理
		foreach my $isoform (map {[unpack($bed12_template, $_)]} @{$locus}) {
			# アイソフォームIDを更新
			$isoform_id++;
			
			# ブロック情報を変換
			$isoform->[10] = [unpack("L*", $isoform->[10])];
			$isoform->[11] = [unpack("L*", $isoform->[11])];
			
			# bed形式の場合
			if ($_[3] eq "bed") {
				$isoform->[3] = "$_[2]$_[1].$isoform_id:" . substr($isoform->[3], index($isoform->[3], ":") + 1);
				$isoform->[5] = $isoform->[5] > 0 ? "+" : "-";
				print join("\t", @{$isoform}[0..9]), "\t", join(",", @{$isoform->[10]}), "\t", join(",", @{$isoform->[11]}), "\n";
			}
			
			# gtf形式の場合
			elsif ($_[3] eq "gtf") {
				# 開始点をgtf形式に修正
				$isoform->[1]++;
				
				# 逆鎖の場合は個々のブロックを逆順に並べ替え
				$isoform->[10] = [reverse(@{$isoform->[10]})] and $isoform->[11] = [reverse(@{$isoform->[11]})] if $isoform->[5] < 0;
				
				# 向きを登録
				my $sign = $isoform->[5];
				$isoform->[5] = $isoform->[5] > 0 ? "+" : "-";
				
				# トランスクリプト行を出力
				print "$isoform->[0]\tfate\ttranscript\t$isoform->[1]\t$isoform->[2]\t$isoform->[4]\t$isoform->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].$isoform_id", '"; transcript_name "', $isoform->[3], '"; transcript_biotype "', $type{$isoform->[8]}, '";', "\n";
				
				# 個々のヒットを示す行を出力
				for (my $i = 0;$i < $isoform->[9];$i++) {
					# エキソン行を出力
					print "$isoform->[0]\tfate\texon\t", $isoform->[1] + $isoform->[11]->[$i], "\t", $isoform->[1] + $isoform->[10]->[$i] + $isoform->[11]->[$i] - 1, "\t.\t$isoform->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].$isoform_id", '"; exon_number "', $i + 1, '"; transcript_name "', $isoform->[3], '"; transcript_biotype "', $type{$isoform->[8]}, '"; exon_id "', "$_[2]$_[1].$isoform_id.", $i + 1, '";', "\n";
					
					# コーディング領域以外を除外
					next if $type{$isoform->[8]} ne "protein_coding";
					
					# コーディング領域行を出力
					print "$isoform->[0]\tfate\tCDS\t", $isoform->[1] + $isoform->[11]->[$i] + 3 * 0 ** ($isoform->[9] - $i + $sign), "\t", $isoform->[1] + $isoform->[10]->[$i] + $isoform->[11]->[$i] - 1 - 3 * 0 ** ($isoform->[9] - $i - $sign), "\t.\t$isoform->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].$isoform_id", '"; exon_number "', $i + 1, '"; transcript_name "', $isoform->[3], '"; transcript_biotype "', $type{$isoform->[8]}, '"; protein_id "', "$_[2]$_[1].$isoform_id", 'p";', "\n";
				}
			}
		}
	}
	return(1);
}

# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
