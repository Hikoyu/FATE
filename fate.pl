#!/usr/bin/env perl

# Copyright (c) 2016 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "fate.pl";	# ソフトウェアの名前
my $version = "ver.1.0.0";	# ソフトウェアのバージョン
my $note = "FATE is Framework for Annotating Translatable Exons.\n  This software annotates protein-coding genes by a classical homology-based method.";	# ソフトウェアの説明
my $usage = "<required items> [optional items]";	# ソフトウェアの使用法 (コマンド非使用ソフトウェアの時に有効)
### 編集範囲 終了 ###

# コマンドを定義
my %command;
### 編集範囲 開始 ###
$command{"search"} = "Search protein-coding genes under specified conditions";
$command{"filter"} = "Filter already annotated protein-coding genes under specified conditions";
### 編集範囲 終了 ###
my @command_list = sort(keys(%command));

# 指定されたコマンドを確認
my $specified_command = "";
if (@command_list and @ARGV) {
	$specified_command = shift(@ARGV);
	if (!grep {$_ eq $specified_command} @command_list) {&exception::error("unknown command: $specified_command");}
}

# 共通オプションを定義
my %option;
### 編集範囲 開始 ###
$option{"b PATH "} = "Path to blast output file (-outfmt 6) to reuse";
$option{"f STR "} = "Output format <bed|gtf> [bed]";
$option{"n STR "} = "Prefix of each locus name [locus]";
$option{"p INT "} = "Number of parallel processes <1-> [1]";
$option{"t INT "} = "Output biotype (+1:functional genes, +2:truncated genes, +4:pseudogenes) <1-7> [7]";
$option{"v INT "} = "Maximum number of candidate isoforms <1->";
### 編集範囲 終了 ###

# コマンドごとのオプション定義を取得
if ($specified_command) {&{\&{"${specified_command}::define"}};}
my @option_list = sort(keys(%option));

# ヘルプを表示 (引数未指定時)
if (!@ARGV and !-p STDIN) {&exception::help;}

# オプションの入力処理
my %opt;
$_ = join("", @option_list);
$_ =~ s/\s+\S+\s+/:/g;
getopts($_, \%opt);

# 未指定オプションのデフォルト値を入力
foreach (@option_list) {
	if ($option{$_} =~ /\[.+\]$/ and !defined($opt{substr($_, 0, 1)})) {
		$opt{substr($_, 0, 1)} = substr($option{$_}, index($option{$_}, "[") + 1, index($option{$_}, "]") - index($option{$_}, "[") - 1);
	}
}

### 編集範囲 開始 ###
# 追加のモジュールを宣言
use List::Util;
use IPC::Open2;
use IO::Pipe;

# 依存するソフトウェアがインストールされているか確認
if (!`which makeblastdb` or !`which blastn` or !`which blastx` or !`which tblastn`) {&exception::error("blast+ not installed");}

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

# 処理を追加
### 編集範囲 終了 ###

# メインルーチンを実行
&main;
exit 0;

# メインルーチン
sub main {
	### 編集範囲 開始 ###
	# 指定された共通オプションを確認
	$opt{"f"} = lc($opt{"f"});
	if ($opt{"f"} ne "bed" and $opt{"f"} ne "gtf") {&exception::error("specify bed or gtf: -f $opt{f}");}
	if ($opt{"p"} !~ /^\d+$/ or $opt{"p"} < 1) {&exception::error("specify INT >= 1: -p $opt{p}");}
	if ($opt{"t"} != 1 and $opt{"t"} != 2 and $opt{"t"} != 3 and $opt{"t"} != 4 and $opt{"t"} != 5 and $opt{"t"} != 6 and $opt{"t"} != 7) {&exception::error("specify INT 1-7: -t $opt{t}");}
	$opt{"t"} = sprintf("%03b", $opt{"t"});
	if (defined($opt{"v"}) and ($opt{"v"} !~ /^\d+$/ or $opt{"v"} < 1)) {&exception::error("specify INT >= 1: -v $opt{v}");}
	if (!defined($opt{"v"})) {$opt{"v"} = 0;}
	
	# 相同性検索の出力ファイルを確認 (-b指定時)
	if ($opt{"b"}) {
		if (!-f $opt{"b"}) {&exception::error("file not found: $opt{b}");}
		if (!-r $opt{"b"}) {&exception::error("file unreadable: $opt{b}");}
		if (!-s $opt{"b"}) {&exception::error("null file specified: $opt{b}");}
	}
	### 編集範囲 終了 ###
	
	# コマンドの実行 (コマンド指定時)
	if ($specified_command) {&{\&{"${specified_command}::body"}};}
	
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
	if ($specified_command) {print STDERR $specified_command;} else {print STDERR $version;}
	print STDERR "\n\nFunctions:\n  $note\n\nUsage:\n  $software ";
	if (!$specified_command and @command_list) {
		print STDERR "<command>\n";
		print STDERR "\nCommand:\n";
		foreach (@command_list) {print STDERR "  $_\t$command{$_}\n";}
	} else {
		if ($specified_command) {print STDERR "$specified_command ";}
		if (@option_list) {print STDERR "[options] ";}
		print STDERR "$usage\n";
		if (@option_list) {
			print STDERR "\nOptions:\n";
			foreach (@option_list) {print STDERR "  -$_\t$option{$_}\n";}
		}
	}
	exit 0;
}

# エラー表示
sub error {
	print STDERR $software;
	if ($specified_command) {print STDERR " $specified_command";}
	print STDERR ": Error: $_[0]";
	if ($_[1]) {print STDERR ": $_[1] line $.";}
	print STDERR "\n";
	exit 1;
}

# 注意表示
sub caution {
	print STDERR $software;
	if ($specified_command) {print STDERR " $specified_command";}
	print STDERR ": Caution: $_[0]";
	if ($_[1]) {print STDERR ": $_[1] line $.";}
	print STDERR "\n";
	return(1);
}

### 編集範囲 開始 ###
## ここからsearchコマンドのパッケージ ##
package search;

# コマンドとオプションを定義
sub define {
	$note = "Search protein-coding genes under specified conditions.";
	$usage = "<genome.fa> <STDIN | in1.fa> [in2.fa ...] [> out.bed | > out.gtf]";
	$option{"s"} = "\tForce GT-AG rule for splice junctions";
	$option{"g STR "} = "Gene prediction program <exonerate|genewise>";
	$option{"h STR "} = "Homology search program <tblastn|blastn|megablast|dc-megablast> [tblastn]";
	$option{"5 INT "} = "Length of 5' flanking region <0-> [300]";
	$option{"3 INT "} = "Length of 3' flanking region <0-> [300]";
	$option{"i INT "} = "Maximum interval length allowed to assemble initial hits [100000]";
	$option{"o INT "} = "Maximum overlap/gap length of query boundries allowed to assemble initial hits [30]";
	$option{"l INT "} = "Minimum length to regard as complete CDS <0-> [0]";
	$option{"c NUM "} = "Minimum query coverage to regard as complete CDS <0-1> [0.85]";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	if ($opt{"g"} and $opt{"g"} ne "exonerate" and $opt{"g"} ne "genewise") {&exception::error("unknown program specified: -g $opt{g}");}
	if ($opt{"h"} ne "tblastn" and $opt{"h"} ne "blastn" and $opt{"h"} ne "megablast" and $opt{"h"} ne "dc-megablast") {&exception::error("unknown program specified: -h $opt{h}");}
	if ($opt{"s"} and !$opt{"g"}) {&exception::caution("-s ignored under -g unspecified");}
	if ($opt{"5"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -5 $opt{5}");}
	if ($opt{"3"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -3 $opt{3}");}
	if ($opt{"i"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -i $opt{i}");}
	if ($opt{"o"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -o $opt{o}");}
	if ($opt{"l"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -l $opt{l}");}
	if ($opt{"c"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"c"} > 1) {&exception::error("specify NUM 0-1: -c $opt{c}");}
	
	# ゲノム配列のfastaファイル名を取得
	my $genome_file = shift(@ARGV);
	
	# 入力ファイルを確認
	if (!@ARGV and !-p STDIN) {&exception::error("input file not specified");}
	foreach (@ARGV) {
		if (!-f $_) {&exception::error("file not found: $_");}
		if (!-r $_) {&exception::error("file unreadable: $_");}
		if (!-s $_) {&exception::error("null file specified: $_");}
	}
	
	# ゲノム配列のfastaファイルを検索
	&common::find_db($genome_file);
	
	# ゲノム配列のfastaインデックスを取得
	my $genome_faidx = &common::read_fasta($genome_file);
	
	# ゲノム配列のblastデータベースを確認 (-b未指定時)
	if (!$opt{"b"}) {&common::check_blastdb($genome_file);}
	
	# 処理を定義
	my $search_engine = $opt{"h"} eq "tblastn" ? "tblastn" : "blastn";
	my $homology_search = "$search_engine -num_threads $opt{p} -outfmt 6 -soft_masking true -db $genome_file";
	$homology_search .= $opt{"h"} eq "tblastn" ? "" : " -task $opt{h}";
	my $gene_prediction = $opt{"g"};
	if ($opt{"g"} and $opt{"g"} eq "exonerate") {
		$gene_prediction .= " -m protein2genome -V 0 --showtargetgff T --showalignment F --showvulgar F";
		$gene_prediction .= $opt{"s"} ? " --forcegtag T" : "";
	}
	elsif ($opt{"g"} and $opt{"g"} eq "genewise") {
		$gene_prediction .= " -divide '#' -pseudo -sum -gff -silent -quiet";
		$gene_prediction .= $opt{"s"} ? "" : " -nosplice_gtag";
	}
	if (`which tee`) {$homology_search .= " | tee fate_search_${search_engine}.out";}
	
	# フォルダを作成 (-g指定時)
	if ($opt{"g"}) {
		if (!-d "query") {mkdir("query") or &exception::error("failed to make directory: query");}
		if (!-d "loci") {mkdir("loci") or &exception::error("failed to make directory: loci");}
	}
	
	# プロセス間通信のパイプを作成
	my $output = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid)) {&exception::error("failed to fork process");}
	
	## ここから子プロセスの処理 ##
	if (!$pid) {
		# 変数を宣言
		my @pid = ();
		my @pipe = ();
		my $pflag = 0;
		
		# プロセス間通信のパイプを作成
		my $order = IO::Pipe->new;
		
		# 指定したプロセス数で並列処理 (-g指定時)
		for (my $pnum = 0;$pnum < $opt{"p"} and $opt{"g"};$pnum++) {
			# プロセス間通信のパイプを作成
			push(@pipe, my $input = IO::Pipe->new);
			
			# プロセス分岐
			$pid[$pnum] = fork;
			
			# プロセス分岐に失敗した場合
			if (!defined($pid[$pnum])) {$pflag = 1;last;}
			
			## ここから孫プロセスの処理 ##
			if (!$pid[$pnum]) {
				# 相同性検索の出力を閉じる
				close(SEARCH_OUT);
				
				# パイプを開く
				$output->writer;
				$order->writer;
				$input->reader;
				
				# ゲノム配列のfastaファイルを開く
				open(GENOME, "<", $genome_file) or &exception::error("failed to open file: $genome_file");
				
				# 子プロセスにデータ要求を送信
				syswrite($order, "$pnum\n");
				
				# 子プロセスから領域データを受信し、遺伝子予測を実行
				while (<$input>) {
					# 改行コードを除去
					chomp;
					
					# タブ文字でデータを分割
					my @col = split(/\t/);
					
					# フランキング配列長に合わせて開始点と終止点を修正
					my $locus_start = $col[5] eq "+" ? $col[1] - $opt{"5"} : $col[1] - $opt{"3"};
					my $locus_end = $col[5] eq "+" ? $col[2] + $opt{"3"} : $col[2] + $opt{"5"};
					if ($locus_start < 0) {$locus_start = 0;}
					if ($locus_end > $genome_faidx->{$col[0]}->{"seq_length"}) {$locus_end = $genome_faidx->{$col[0]}->{"seq_length"};}
					
					# 取り出す領域の位置情報を算出
					my $root = $genome_faidx->{$col[0]};
					my $start_point = int($locus_start / $root->{"row_width"}) * $root->{"row_bytes"} + $locus_start % $root->{"row_width"};
					my $end_point = int($locus_end / $root->{"row_width"}) * $root->{"row_bytes"} + $locus_end % $root->{"row_width"};
					
					# 領域の塩基配列を取得
					seek(GENOME, $root->{"seq_start"} + $start_point, 0);
					read(GENOME, my $locus_seq, $end_point - $start_point);
					
					# 配列に含まれる改行コードを除去
					$locus_seq =~ s/\n|\r//g;
					
					# 配列を大文字に変換
					$locus_seq = uc($locus_seq);
					
					# 相補鎖に変換
					if ($col[5] eq "-") {&common::complementary($locus_seq);}
					
					# 領域開始点を表示形式に合わせる
					$locus_start++;
					
					# 領域の配列ファイルを作成
					open(LOCUS, ">", "loci/locus$pnum.fa") or &exception::error("failed to make file: loci/locus$pnum.fa");
					
					# 取得した配列をファイルに出力
					print LOCUS ">$col[0]:$locus_start-$locus_end($col[5])\n";
					for (my $pos = 0;$pos < length($locus_seq);$pos += 60) {print LOCUS substr($locus_seq, $pos, 60), "\n";}
					
					# 領域の配列ファイルを閉じる
					close(LOCUS);
					
					# 変数を宣言
					my @genes = ();
					my $query_start = 0;
					my $query_end = 0;
					my $summary_flag = 0;
					
					# 遺伝子予測を実行
					my $args = $opt{"g"} eq "exonerate" ? "-q query/$col[3].fa -t loci/locus$pnum.fa" : $opt{"g"} eq "genewise" ? "query/$col[3].fa loci/locus$pnum.fa" : "";
					open(PREDICT, "-|", "$gene_prediction $args 2>/dev/null") or &exception::error("failed to execute gene prediction: $col[0]:$locus_start-$locus_end($col[5]) vs $col[3]");
					
					# 領域開始点をbed形式に戻す
					$locus_start--;
					
					# 遺伝子予測結果を1行ずつ読み込んで処理
					while (<PREDICT>) {
						# コメント行を除外
						if (substr($_, 0, 1) eq "#") {$summary_flag = 1;next;}
						
						# summary行の処理 (-g genewise指定時)
						if (!$summary_flag and $opt{"g"} eq "genewise") {
							# 空白文字でデータを分割
							my @col = split(/\s+/);
							if ($col[0] eq "Bits") {next;}
							$query_start = ($col[2] - 1) * 3;
							$query_end = $col[3] * 3;
							next;
						}
						
						# タブ文字でデータを分割
						my @gff = split(/\t/);
						
						# 逆向きの予測を除外
						if ($gff[6] eq "-") {next;}
						
						# 開始点をbed形式に合わせる
						$gff[3]--;
						
						# match行またはgene行の処理
						if ($gff[2] eq "match" or $gff[2] eq "gene") {push(@genes, [$col[0], $gff[3], $gff[4], $col[3], $gff[5], "+", 0, 0, "", 0, [], [], $query_start, $query_end]);}
						
						# cds行の処理
						elsif ($gff[2] eq "cds") {$genes[-1]->[9]++;push(@{$genes[-1]->[10]}, $gff[4] - $gff[3]);push(@{$genes[-1]->[11]}, $gff[3]);}
						
						# similarity行の処理 (-g exonerate指定時)
						elsif ($gff[2] eq "similarity") {
							my @query_block = ();
							my @col = split(/ ; /, $gff[8]);
							foreach (@col) {
								my ($key, @summary) = split(/ /);
								if ($key eq "Align") {push(@query_block, (($summary[1] - 1) * 3, ($summary[1] - 1) * 3 + $summary[2]));}
							}
							$genes[-1]->[12] = List::Util::min(@query_block);
							$genes[-1]->[13] = List::Util::max(@query_block);
						}
					}
					
					# 遺伝子予測を終了
					close(PREDICT);
					
					# 翻訳領域を推定
					foreach my $gene (@genes) {
						# 変数を宣言
						my $upstream_truncation = 0;
						my $downstream_truncation = 0;
						my $completeness = 0;
						
						# Nを含まない上流配列を取得し、その長さが指定値より短い場合は5'側truncateとする
						my $upstream_seq = substr(substr($locus_seq, 0, $gene->[1]), -$opt{"5"});
						my $upstream_offset = $col[5] eq "+" ? $col[1] - $locus_start - $gene->[1] : $locus_end - $gene->[1] - $col[2];
						if ($upstream_offset < 0) {$upstream_offset = 0;}
						$upstream_seq =~ s/^.*N//;
						if (length($upstream_seq) < $opt{"5"} - $upstream_offset) {$upstream_truncation = 1;}
						
						# Nを含まない下流配列を取得し、その長さが指定値より短い場合は3'側truncateとする
						my $downstream_seq = substr($locus_seq, $gene->[2], $opt{"3"});
						my $downstream_offset = $col[5] eq "+" ? $locus_start + $gene->[2] - $col[2] : $col[1] - $locus_end + $gene->[2];
						if ($downstream_offset < 0) {$downstream_offset = 0;}
						$downstream_seq =~ s/N.*$//;
						if (length($downstream_seq) < $opt{"3"} - $downstream_offset) {$downstream_truncation = 1;}
						
						# 5'側クエリー被覆率が指定値以上の場合
						if ($col[4] - $gene->[12] >= $col[4] * $opt{"c"}) {
							# 上流フランキング領域のアミノ酸配列を取得
							my $upstream_aa = &common::translate($upstream_seq, length($upstream_seq) % 3);
							
							# 先頭エキソンのアミノ酸配列を取得
							my $first_exon_aa = &common::translate(substr($locus_seq, $gene->[11]->[0], $gene->[10]->[0]), 0);
							
							# 開始コドンを探索
							my $start_pos = rindex($upstream_aa, "M");
							if (rindex($upstream_aa, "*") < $start_pos) {$start_pos -= length($upstream_aa);}
							else {$start_pos = index($first_exon_aa, "M");$start_pos += 0 ** ($start_pos + 1);}
							
							# 開始コドンを考慮しても5'側クエリー被覆率が指定値以上の場合は先頭ブロックを修正
							if ($col[4] - $gene->[12] - $start_pos * 3 >= $col[4] * $opt{"c"}) {$gene->[10]->[0] -= $start_pos * 3;$gene->[11]->[0] += $start_pos * 3;}
							
							# 5'端が完全であることを認定
							$upstream_truncation = 0;
							$completeness++;
						}
						
						# 3'側クエリー被覆率が指定値以上の場合
						if ($gene->[13] >= $col[4] * $opt{"c"}) {
							# 下流フランキング領域のアミノ酸配列を取得
							my $downstream_aa = &common::translate($downstream_seq, 0);
							
							# 末尾エキソンのアミノ酸配列を取得
							my $last_exon_aa = &common::translate(substr($locus_seq, $gene->[11]->[-1], $gene->[10]->[-1]), $gene->[10]->[-1] % 3);
							
							# 終止コドンを探索
							my $terminal_pos = rindex($last_exon_aa, "*");
							if ($terminal_pos >= 0) {$terminal_pos -= length($last_exon_aa);}
							else {$terminal_pos = index($downstream_aa, "*");}
							$terminal_pos++;
							
							# 終止コドンを考慮しても3'側クエリー被覆率が指定値以上の場合は末尾ブロックを修正
							if ($gene->[13] + $terminal_pos * 3 >= $col[4] * $opt{"c"}) {$gene->[10]->[-1] += $terminal_pos * 3;}
							
							# 3'端が完全であることを認定
							$downstream_truncation = 0;
							$completeness++;
						}
						
						# 領域開始点・終止点を修正
						($gene->[1], $gene->[2]) = ($gene->[11]->[0], $gene->[10]->[-1] + $gene->[11]->[-1]);
						($gene->[6], $gene->[7]) = ($gene->[1], $gene->[2]);
						
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
						elsif ($completeness < 2) {$gene->[8] = $upstream_truncation + $downstream_truncation ? "yellow" : "red";}
						
						# 配列長が指定値未満の場合 (偽遺伝子)
						elsif (List::Util::sum(@{$gene->[10]}) < $opt{"l"}) {$gene->[8] = "red";}
						
						# 上記に該当しない場合 (機能遺伝子)
						else {$gene->[8] = "blue";}
						
						# ゲノム配列の座標で修正
						($gene->[1], $gene->[2]) = ($col[5] eq "+" ? $locus_start + $gene->[1] : $locus_end - $gene->[2], $col[5] eq "+" ? $locus_start + $gene->[2] : $locus_end - $gene->[1]);
						($gene->[6], $gene->[7]) = ($gene->[1], $gene->[2]);
						$gene->[5] = $col[5];
						if ($col[5] eq "-") {
							$gene->[10] = [reverse(@{$gene->[10]})];
							$gene->[11] = [reverse(@{$gene->[11]})];
							for (my $i = 0;$i < @{$gene->[11]};$i++) {$gene->[11]->[$i] = $gene->[2] - $gene->[1] - $gene->[10]->[$i] - $gene->[11]->[$i];}
						}
						
						# bed12形式で親プロセスに送信
						syswrite($output, join("\t", (@{$gene}[0..9], join(",", @{$gene->[10]}), join(",", @{$gene->[11]}))) . "\n");
					}
					
					# 親プロセスに染色体、連鎖群、スキャフォールドまたはコンティグ名を送信
					syswrite($output, "$col[0]\n");
					
					# 子プロセスにプロセス番号を送信
					syswrite($order, "$pnum\n");
				}
				
				# ゲノム配列のfastaファイルを閉じる
				close(GENOME);
				
				# プロセスを終了
				exit 0;
			}
			## ここまで孫プロセスの処理 ##
		}
		
		# プロセス分岐に失敗していた場合
		if ($pflag) {&exception::error("failed to fork process");}
		
		# 相同性検索の出力ファイルを開く (-b指定時)
		if ($opt{"b"}) {open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}");}
		
		# 相同性検索を実行
		else {main::open2(*SEARCH_OUT, *SEARCH_IN, $homology_search) or &exception::error("failed to execute homology search");}
		
		# 変数を宣言
		my %query_len = ();
		
		# プロセス間通信のパイプを作成
		my $query_report = IO::Pipe->new;
		
		# プロセス分岐
		my $pid = fork;
		
		# プロセス分岐に失敗した場合
		if (!defined($pid)) {&exception::error("failed to fork process");}
		
		## ここから孫プロセスの処理 ##
		if (!$pid) {
			# 相同性検索の出力を閉じる
			close(SEARCH_OUT);
			
			# 変数を宣言
			my $query_title = "";
			my $query_seq = "";
			
			# クエリー配列を読み込みながら処理
			while (my $line = <>) {
				# 改行コードを除去
				$line =~ s/\n|\r//g;
				
				# 配列行の処理
				if ($line !~ /^>/) {$query_seq .= uc($line);}
				
				# ID行の処理
				if ($line =~ /^>/ or eof) {
					if ($query_title) {
						# クエリー配列を相同性検索に入力 (-b未指定時)
						if (!$opt{"b"}) {
							print SEARCH_IN ">$query_title\n";
							for (my $pos = 0;$pos < length($query_seq);$pos += 60) {print SEARCH_IN substr($query_seq, $pos, 60), "\n";}
						}
						
						# クエリー配列長をハッシュに登録
						$query_len{$query_title} = length($query_seq);
						
						# クエリー配列をアミノ酸配列に翻訳 (-h tblastn非指定時)
						if ($opt{"h"} ne "tblastn") {$query_seq = &common::translate($query_seq, 0);}
						
						# クエリー配列からアスタリスクを除去
						$query_seq =~ s/\*//g;
						
						# クエリー配列をファイルに出力 (-g指定時)
						if ($opt{"g"}) {
							# 出力ファイルを作成
							open(QUERY, ">", "query/$query_title.fa") or &exception::error("failed to make file: query/$query_title.fa");
							
							# クエリー配列をfasta形式で個別のファイルに出力
							print QUERY ">$query_title\n";
							for (my $pos = 0;$pos < length($query_seq);$pos += 60) {print QUERY substr($query_seq, $pos, 60), "\n";}
							
							# 出力ファイルを閉じる
							close(QUERY);
						}
					}
					
					# ID行の最初の空白文字の前までをタイトルとして登録
					($query_title) = split(/\s/, substr($line, 1));
					
					# 配列をリセット
					$query_seq = "";
				}
			}
			
			# 相同性検索への入力を閉じる (-b未指定時)
			if (!$opt{"b"}) {close(SEARCH_IN);}
			
			# パイプを開く
			$query_report->writer;
			
			# 子プロセスにクエリー配列長データを送信
			foreach (keys(%query_len)) {syswrite($query_report, "$_\t$query_len{$_}\n");}
			
			# プロセスを終了
			exit 0;
		}
		## ここまで孫プロセスの処理 ##
		
		# 相同性検索への入力を閉じる (-b未指定時)
		if (!$opt{"b"}) {close(SEARCH_IN);}
		
		# 変数を宣言
		my %blast_hits = ();
		my $num_hits = 0;
		
		# 相同性検索の出力を読み込みながら処理
		print STDERR "Running homology search...";
		while (<SEARCH_OUT>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# データを保存
			my $dat = $col[8] < $col[9] ? {"query_start" => $col[6] - 1, "query_end" => $col[7], "locus_start" => $col[8] - 1, "locus_end" => $col[9], "strand" => 1}
										: {"query_start" => $col[7], "query_end" => $col[6] - 1, "locus_start" => $col[9] - 1, "locus_end" => $col[8], "strand" => -1};
			$dat->{"score"} = $col[11];
			$dat->{"total_score"} = $col[11];
			$dat->{"locus_destination"} = $dat->{"locus_end"};
			$dat->{"block_size"} = [$dat->{"locus_end"} - $dat->{"locus_start"}];
			$dat->{"block_start"} = [$dat->{"locus_start"}];
			$dat->{"num_connection"} = 0;
			push(@{$blast_hits{$col[1]}{$col[0]}}, $dat);
			
			# ヒット数を加算
			$num_hits++;
		}
		print STDERR "completed\n";
		
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# パイプを開く
		$query_report->reader;
		
		# 孫プロセスからクエリー配列長データを受信
		print STDERR "Scanning queries...";
		while (<$query_report>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# クエリー配列長をハッシュに登録
			$query_len{$col[0]} = $col[1];
		}
		print STDERR "completed\n";
		
		# 孫プロセスを刈り取る
		waitpid($pid, 0);
		if ($?) {&exception::error("process abnormally exited");}
		
		# 変数を宣言
		my $fin_hits = 0;
		
		# パイプを開く
		if ($opt{"g"}) {$order->reader;}
		for (my $pnum = 0;$pnum < @pipe;$pnum++) {$pipe[$pnum]->writer;}
		$output->writer;
		
		# 相同性検索の各ヒット領域を整理し、領域データを送信
		print STDERR "Running gene prediction...";
		foreach my $subject (sort(keys(%blast_hits))) {
			# 変数を宣言
			my $num_loci = 0;
			
			foreach my $query (keys(%{$blast_hits{$subject}})) {
				# ヒットを領域開始点と領域終止点の順で並べ替え
				my @sorted_hits = sort {$a->{"locus_start"} <=> $b->{"locus_start"} or $a->{"locus_end"} <=> $b->{"locus_end"}} @{$blast_hits{$subject}{$query}};
				
				# 後方のヒットから判定
				for (my $i = -1;$i >= -@sorted_hits;$i--) {
					# 連結するヒットを決定
					my $k = $i + 1;
					while ($k < 0) {
						# 指定値を超える距離のヒットに到達した時点で終了
						if ($sorted_hits[$k]->{"locus_start"} - $sorted_hits[$i]->{"locus_end"} > $opt{"i"}) {last;}
						
						# 方向が異なるヒットを除外
						if (!($sorted_hits[$k]->{"strand"} + $sorted_hits[$i]->{"strand"})) {next;}
						
						# 指定値を超えるクエリーオーバーラップのヒットを除外
						if (abs($sorted_hits[$k]->{"query_start"} - $sorted_hits[$i]->{"query_end"}) > $opt{"o"}) {next;}
						
						# クエリー開始点が前進しないヒットを除外
						if ($sorted_hits[$k]->{"query_start"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_start"} * $sorted_hits[$i]->{"strand"}) {next;}
						
						# クエリー終止点が前進しないヒットを除外
						if ($sorted_hits[$k]->{"query_end"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_end"} * $sorted_hits[$i]->{"strand"}) {next;}
						
						# 総スコアが更新されないヒットを除外
						if ($sorted_hits[$k]->{"total_score"} <= $sorted_hits[$i]->{"total_score"} - $sorted_hits[$i]->{"score"}) {next;}
						
						# データを更新
						$sorted_hits[$i]->{"total_score"} = $sorted_hits[$k]->{"total_score"} + $sorted_hits[$i]->{"score"};
						$sorted_hits[$i]->{"locus_destination"} = $sorted_hits[$k]->{"locus_destination"};
						$sorted_hits[$i]->{"block_size"} = [$sorted_hits[$i]->{"locus_end"} - $sorted_hits[$i]->{"locus_start"}, @{$sorted_hits[$k]->{"block_size"}}];
						$sorted_hits[$i]->{"block_start"} = [$sorted_hits[$i]->{"locus_start"}, @{$sorted_hits[$k]->{"block_start"}}];
						$sorted_hits[$k]->{"num_connection"}++;
					}
					continue {
						# 領域の区切りを判定
						if ($sorted_hits[$k]->{"strand"} + $sorted_hits[$i]->{"strand"} and !$sorted_hits[$k]->{"num_connection"}) {last;}
						$k++;
					}
				}
				
				# 孫プロセスにデータを送信 (-g指定時)
				if ($opt{"g"}) {
					foreach my $locus (@sorted_hits) {
						# 開始ブロックでないヒットを除外
						if ($locus->{"num_connection"}) {next;}
						
						# プロセス番号を受信した孫プロセスにbed6形式でデータを送信
						syswrite($pipe[<$order>], join("\t", ($subject, $locus->{"locus_start"}, $locus->{"locus_destination"}, $query, $opt{"h"} eq "tblastn" ? $query_len{$query} * 3 : $query_len{$query}, $locus->{"strand"} > 0 ? "+" : "-")) . "\n");
						
						# 領域数を加算
						$num_loci++;
					}
				}
				
				# 親プロセスにデータを送信 (-g未指定時)
				else {
					foreach my $locus (@sorted_hits) {
						# 開始ブロックでないヒットを除外
						if ($locus->{"num_connection"}) {next;}
						
						# 親プロセスにbed12形式でデータを送信
						syswrite($output, join("\t", ($subject, $locus->{"locus_start"}, $locus->{"locus_destination"}, $query, $locus->{"score"}, $locus->{"strand"} > 0 ? "+" : "-", $locus->{"locus_start"}, $locus->{"locus_destination"}, ".", scalar(@{$locus->{"block_size"}}), join(",", @{$locus->{"block_size"}}), join(",", map {$_ - $locus->{"locus_start"}} @{$locus->{"block_start"}}))) . "\n$subject\n");
						
						# 領域数を加算
						$num_loci++;
					}
				}
				
				# 送信済みヒット数を加算
				$fin_hits += @sorted_hits;
				
				# 途中経過を表示
				print STDERR "\rRunning gene prediction...",int($fin_hits / $num_hits * 100),"%";
			}
			
			# 親プロセスに染色体、連鎖群、スキャフォールドまたはコンティグ名と領域数を送信
			syswrite($output, "$subject\t$num_loci\n");
			
			# 処理が完了したデータを削除
			delete($blast_hits{$subject});
		}
		
		# パイプを閉じる
		undef(@pipe);
		
		# 孫プロセスを刈り取る
		foreach (@pid) {
			waitpid($_, 0);
			if ($?) {&exception::error("process abnormally exited");}
		}
		print STDERR "\rRunning gene prediction...completed\n";
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスの処理 ##
	
	# 変数を宣言
	my %loci = ();
	my %num_loci = ();
	my $locus_num = 0;
	
	# パイプを開く
	$output->reader;
	
	# データを受信しながら処理
	while (<$output>) {
		# 改行コードを除去
		chomp;
		
		# タブ文字でデータを分割
		my @col = split(/\t/);
		
		# bed形式のデータを受信した場合
		if (@col > 2) {push(@{$loci{$col[0]}}, \@col);next;}
		
		# 染色体、連鎖群、スキャフォールドまたはコンティグ名と領域数を受信した場合
		if (@col > 1) {$num_loci{$col[0]}->[0] = $col[1];}
		
		# 染色体、連鎖群、スキャフォールドまたはコンティグ名のみ受信した場合
		else {$num_loci{$col[0]}->[1]++;}
		
		# 各染色体、連鎖群、スキャフォールドまたはコンティグについて順に処理
		foreach my $subject (sort(keys(%num_loci))) {
			# 領域数が一致しない場合は以降の処理も含めて保留
			if (!$num_loci{$subject}->[0] or !$num_loci{$subject}->[1] or $num_loci{$subject}->[0] > $num_loci{$subject}->[1]) {last;}
			
			# データを出力
			&common::modify_bed($loci{$subject}, $locus_num, $opt{"n"}, $opt{"t"}, $opt{"v"}, $opt{"f"});
			
			# 処理が完了したデータを削除
			delete($loci{$subject});
			delete($num_loci{$subject});
		}
	}
	
	# 子プロセスを刈り取る
	waitpid($pid, 0);
	if ($?) {&exception::error("process abnormally exited");}
	
	return(1);
}

## ここからfilterコマンドのパッケージ ##
package filter;

# コマンドとオプションを定義
sub define {
	$note = "Filter already annotated protein-coding genes under specified conditions.";
	$usage = "<STDIN | in1.bed> [in2.bed ...] [> out.bed | > out.gtf]";
	$option{"d PATH "} = "Path to gene or protein database file (fasta format)";
	$option{"g PATH "} = "Path to target genome data file (fasta format)";
	$option{"h STR "} = "Homology search program <blastx|blastn|megablast|dc-megablast> [blastx]";
	$option{"k STR "} = "Keywords for filtering (AND[&], OR[;], BUT[!])";
	$option{"r INT "} = "Cutoff rank of hits <1->";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	if ($opt{"b"} and $opt{"d"}) {&exception::error("options incompatible: -b and -d");}
	if (!$opt{"g"} and $opt{"d"}) {&exception::error("-g required under -d");}
	if ($opt{"h"} ne "blastx" and $opt{"h"} ne "blastn" and $opt{"h"} ne "megablast" and $opt{"h"} ne "dc-megablast") {&exception::error("unknown program specified: -h $opt{h}");}
	if (defined($opt{"r"}) and ($opt{"r"} !~ /^\d+$/ or $opt{"r"} < 1)) {&exception::error("specify INT >= 1: -r $opt{r}");}
	
	# 入力ファイルを確認
	if (!@ARGV and !-p STDIN) {&exception::error("input file not specified");}
	foreach (@ARGV) {
		if (!-f $_) {&exception::error("file not found: $_");}
		if (!-r $_) {&exception::error("file unreadable: $_");}
		if (!-s $_) {&exception::error("null file specified: $_");}
	}
	
	# 処理を定義
	my $search_engine = $opt{"h"} eq "blastx" ? "blastx" : "blastn";
	my $homology_search = "$search_engine -num_threads $opt{p} -outfmt '6 std salltitles' -soft_masking true -strand plus";
	$homology_search .= $opt{"h"} eq "blastx" ? "" : " -task $opt{h}";
	if ($opt{"d"}) {$homology_search .= " -db $opt{d}";}
	if (`which tee`) {$homology_search .= " | tee fate_filter_${search_engine}.out";}

	# 変数を宣言
	my $genome_faidx = "";
	
	# ゲノム配列のfastaファイルを確認 (-g指定時)
	if ($opt{"g"}) {
		# ゲノム配列のfastaファイルを検索
		&common::find_db($opt{"g"});
		
		# ゲノム配列のfastaインデックスを取得
		$genome_faidx = &common::read_fasta($opt{"g"});
	}
	
	# 参照データファイルを確認 (-d指定時)
	if ($opt{"d"}) {
		# 参照配列のfastaファイルを検索
		&common::find_db($opt{"d"});
		
		# 参照配列のblastデータベースを確認 (-b未指定時)
		if ($opt{"h"} eq "blastx") {&common::check_blastdb($opt{"d"}, 1);} else {&common::check_blastdb($opt{"d"});}
	}
	
	# プロセス間通信のパイプを作成
	my $output = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid)) {&exception::error("failed to fork process");}
	
	## ここから子プロセスの処理 ##
	if (!$pid) {
		# パイプを開く
		$output->writer;
		
		# 相同性検索を実行 (-d指定時)
		if ($opt{"d"}) {main::open2(*SEARCH_OUT, *SEARCH_IN, $homology_search) or &exception::error("failed to execute homology search");}
		
		# 相同性検索の出力ファイルを開く (-b指定時)
		elsif ($opt{"b"}) {open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}");}
		
		# 相同性検索を実行せずに処理
		else {
			# 入力データを親プロセスに送信
			while (<>) {syswrite($output, $_);}
			
			# プロセスを終了
			exit 0;
		}
		
		# プロセス間通信のパイプを作成
		my $input = IO::Pipe->new;
		
		# プロセス分岐
		my $pid = fork;
		
		# プロセス分岐に失敗した場合
		if (!defined($pid)) {&exception::error("failed to fork process");}
		
		## ここから孫プロセスの処理 ##
		if (!$pid) {
			# 相同性検索の出力を閉じる
			close(SEARCH_OUT);
			
			# 変数を宣言
			my @bed = ();
			
			# 入力データをリストに保存
			while (<>) {push(@bed, $_);}
			
			# 配列に変換して相同性検索に入力 (-d指定時)
			if ($opt{"d"}) {
				# ゲノム配列のfastaファイルを開く
				open(GENOME, "<", $opt{"g"}) or &exception::error("failed to open file: $opt{g}");
				
				# 1領域ずつ処理
				foreach (@bed) {
					# タブ文字でデータを分割
					my @col = split(/\t/);
					
					# 取り出す領域の位置情報を算出
					my $root = $genome_faidx->{$col[0]};
					my $start_point = int($col[1] / $root->{"row_width"}) * $root->{"row_bytes"} + $col[1] % $root->{"row_width"};
					my $end_point = int($col[2] / $root->{"row_width"}) * $root->{"row_bytes"} + $col[2] % $root->{"row_width"};
					
					# 領域の塩基配列を取得
					seek(GENOME, $root->{"seq_start"} + $start_point, 0);
					read(GENOME, my $locus_seq, $end_point - $start_point);
					
					# 改行コードを除去
					$locus_seq =~ s/\n|\r//g;
					
					# 変数を宣言
					my $query_seq = "";
					
					# ブロックの領域だけを連結
					my @block_length = split(/,/, $col[10]);
					my @block_pos = split(/,/, $col[11]);
					for (my $i = 0;$i < $col[9];$i++) {$query_seq .= substr($locus_seq, $block_pos[$i], $block_length[$i]);}
					
					# 相補鎖に変換
					if ($col[5] eq "-") {&common::complementary($query_seq);}
					
					# 配列を相同性検索に入力
					print SEARCH_IN ">$col[3]\n";
					for (my $pos = 0;$pos < length($query_seq);$pos += 60) {print SEARCH_IN substr($query_seq, $pos, 60), "\n";}
				}
				
				# ゲノム配列のfastaファイルを閉じる
				close(GENOME);
				
				# 相同性検索への入力を閉じる
				close(SEARCH_IN);
			}
			
			# パイプを開く
			$input->writer;
			
			# 読み込んだデータを子プロセスに送信
			while (@bed) {syswrite($input, shift(@bed));}
			
			# プロセスを終了
			exit 0;
		}
		## ここまで孫プロセスの処理 ##
		
		# 相同性検索への入力を閉じる (-d指定時)
		if ($opt{"d"}) {close(SEARCH_IN);}
		
		# 変数を宣言
		my %blast_hits = ();
		my %description = ();
		
		# 相同性検索の出力を読み込みながら処理
		print STDERR "Running homology search...";
		while (<SEARCH_OUT>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# ヒット情報を集計
			push(@{$blast_hits{$col[0]}{$col[1]}}, {"query_start" => $col[6] - 1, "query_end" => $col[7], "locus_start" => $col[8] - 1, "locus_end" => $col[9], "score" => $col[11]});
			$description{$col[1]} = $col[12];
		}
		print STDERR "completed\n";
		
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# パイプを開く
		$input->reader;
		
		# 孫プロセスからデータを受信し、フィルタリングを行う
		print STDERR "Filtering data...";
		while (<$input>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# 変数を宣言
			my %score = ();
			
			# 各ヒットについて処理
			foreach my $target (keys(%{$blast_hits{$col[3]}})) {
				# 変数を宣言
				my @query_array = ();
				my @target_array = ();
				my $total_score = 0;
				my $total_length = 0;
				
				# スコアを算出
				foreach my $prefix (@{$blast_hits{$col[3]}{$target}}) {
					$total_score += $prefix->{"score"};
					$total_length += $prefix->{"query_end"} - $prefix->{"query_start"} + $prefix->{"locus_end"} - $prefix->{"locus_start"};
					@query_array[$prefix->{"query_start"}..$prefix->{"query_end"} - 1] = (1) x ($prefix->{"query_end"} - $prefix->{"query_start"});
					@target_array[$prefix->{"locus_start"}..$prefix->{"locus_end"} - 1] = (1) x ($prefix->{"locus_end"} - $prefix->{"locus_start"});
				}
				$score{$target} = $total_score / $total_length * scalar(grep {$_} (@query_array, @target_array));
			}
			
			# ヒットをスコア順に並べ替え
			my @targets = sort {$score{$b} <=> $score{$a}} keys(%score);

			# 指定値より下位のヒットを削除
			if ($opt{"r"}) {splice(@targets, $opt{"r"});}
			
			# 条件を満たすものを親プロセスに送信
			if (grep {&common::keyword_search($description{$_}, $opt{"k"})} @targets) {syswrite($output, "$_\n");}
		}
		print STDERR "completed\n";
		
		# 孫プロセスを刈り取る
		waitpid($pid, 0);
		if ($?) {&exception::error("process abnormally exited");}
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスの処理 ##
	
	# 変数を宣言
	my @buffer = ();
	my $locus_num = 0;
	my $last_locus = "";
	
	# パイプを開く
	$output->reader;
	
	# 子プロセスからデータを受信しながら処理
	while (<$output>) {
		# 改行コードを除去
		chomp;
		
		# タブ文字でデータを分割
		my @col = split(/\t/);
		
		# 前回の領域と異なる場合
		if ($col[0] ne $last_locus) {
			# 結果を出力
			if ($last_locus) {&common::modify_bed(\@buffer, $locus_num, $opt{"n"}, $opt{"t"}, $opt{"v"}, $opt{"f"});}
			
			# 前回の領域を更新
			$last_locus = $col[0];
			
			# リストを初期化
			@buffer = ();
		}
		
		# 結果をリストに保存
		push(@buffer, \@col);
	}
	
	# 残りの結果を出力
	&common::modify_bed(\@buffer, $locus_num, $opt{"n"}, $opt{"t"}, $opt{"v"}, $opt{"f"});
	
	# 子プロセスを刈り取る
	waitpid($pid, 0);
	if ($?) {&exception::error("process abnormally exited");}
	
	return(1);
}

## ここから共通処理のパッケージ ##
package common;

# キーワード検索 common::keyword_search(検索対象文字列, キーワード)
sub keyword_search {
	# キーワード未指定時は真を返す
	if (!$_[1]) {return(1);}
	
	# 検索を実行して結果を返す
	return(grep {!grep {if (/^!/) {s/^!//;$_[0] =~ /$_/i} else {$_[0] !~ /$_/i}} split(/&/, $_)} split(/;/, $_[1]));
}

# 相補鎖変換 common::complementary(配列)
sub complementary {
	$_[0] = reverse($_[0]);
	$_[0] =~ tr/ATGCRYKMDBVH/TACGYRMKHVBD/;
	return(1);
}

# アミノ酸変換 common::translate(配列, 読み枠)
sub translate {
	my $aa = "";
	for (my $site = $_[1];$site < length($_[0]);$site += 3) {
		my $triplet = substr($_[0], $site, 3);
		if (exists($codon{$triplet})) {$aa .= $codon{$triplet};} else {$aa .= "X";}
	}
	return($aa);
}

# データベース検索 common::find_db(ファイル名)
sub find_db {
	# 環境変数からデータベースパスを取得
	my @db_list = $ENV{"BLASTDB"} ? map {"$_/"} split(/:/, $ENV{"BLASTDB"}) : ();
	
	# ファイルを探索
	my $prefix = List::Util::first {-f "$_$_[0]"} ("", @db_list);
	
	# ファイルが存在しない場合
	if (!defined($prefix)) {&exception::error("file not found: $_[0]");}
	
	# ファイルパスを変更
	$_[0] = "$prefix$_[0]";
	
	# ファイルを確認
	if (!-r $_[0]) {&exception::error("file unreadable: $_[0]");}
	if (!-s $_[0]) {&exception::error("null file specified: $_[0]");}
	return(1);
}

# BLASTデータベース作成 common::check_blastdb(ファイルパス, 配列型)
sub check_blastdb {
	# 配列型を確認
	my $type = "nucl";
	if ($_[1]) {$type = "prot";}
	my $init = substr($type, 0, 1);
	
	# BLASTデータベースファイルを確認
	if ((!-f "$_[0].${init}hr" or !-f "$_[0].${init}in" or !-f "$_[0].${init}sq") and !-f "$_[0].${init}al") {&exception::caution("BLAST database file not found: $_[0].${init}*");}
	elsif ((!-r "$_[0].${init}hr" or !-r "$_[0].${init}in" or !-r "$_[0].${init}sq") and !-r "$_[0].${init}al") {&exception::caution("unreadable BLAST database file: $_[0].${init}*");}
	elsif ((!-s "$_[0].${init}hr" or !-s "$_[0].${init}in" or !-s "$_[0].${init}sq") and !-s "$_[0].${init}al") {&exception::caution("null BLAST database file : $_[0].${init}*");}
	else {return(1);}
	
	# BLASTデータベースファイルを作成
	print STDERR "Building BLAST database...";
	my $check = system("makeblastdb -in $_[0] -dbtype $type 1>/dev/null 2>&1");
	if ($check != 0) {&exception::error("failed to build BLAST database");}
	print STDERR "completed\n";
	return(1);
}

# fastaインデックス読み込み common::read_fasta(ファイルパス)
sub read_fasta {
	# 変数を宣言
	my %fasta_index = ();
	my $id_key = "";
	my $id_start = 0;
	my $seq_length = 0;
	my $seq_start = 0;
	my $row_width = 0;
	my $row_bytes = 0;
	
	# fastaインデックスファイルを確認
	if (!-f "$_[0].fai") {&exception::caution("index file not found: $_[0].fai");}
	elsif (!-r "$_[0].fai") {&exception::caution("unreadable index file: $_[0].fai");}
	elsif (!-s "$_[0].fai") {&exception::caution("null index file: $_[0].fai");}
	
	## ここからインデックスファイル存在時の処理 ##
	else {
		# fastaインデックスファイルを開く
		open(FASTA_INDEX, "<", "$_[0].fai") or &exception::error("failed to open index file: $_[0].fai");
		
		# データを読み込みながら処理
		print STDERR "Loading fasta index...";
		while (my $line = <FASTA_INDEX>) {
			# 改行コードを除去
			$line =~ s/\n|\r//g;
			
			# タブ文字でデータを分割
			($id_key, $seq_length, $seq_start, $row_width, $row_bytes) = split(/\t/, $line);
			
			# データをインデックスハッシュに登録
			$fasta_index{$id_key} = {
				"id_start" => $id_start,
				"seq_length" => $seq_length,
				"seq_start" => $seq_start,
				"row_width" => $row_width,
				"row_bytes" => $row_bytes
			};
			
			# 次のID開始点を登録
			$id_start = $seq_start + int($seq_length / $row_width) * $row_bytes + $seq_length % $row_width;
		}
		print STDERR "completed\n";
		
		# fastaインデックスファイルを閉じる
		close(FASTA_INDEX);
		
		# インデックスハッシュのリファレンスを返す
		return(\%fasta_index);
	}
	
	## ここからインデックスファイル非存在時の処理 ##
	# 変数を宣言
	my $error_flag = 0;
	
	# fastaファイルを開く
	open(FASTA, "<", $_[0]) or &exception::error("failed to open file: $_[0]");
	
	# fastaインデックスファイルを新規作成
	open(FASTA_INDEX, ">", "$_[0].fai") or &exception::error("failed to build index file: $_[0].fai");
	
	# データを読み込みながら処理
	print STDERR "Building fasta index...";
	while (my $line = <FASTA>) {
		# 1行の文字列数を取得
		my $row_bytes1 = length($line);
		
		# 改行コードを除去
		$line =~ s/\n|\r//g;
		
		# 配列行の処理
		if ($id_key and $line !~ /^>/) {
			# データの整合性を確認
			if ($error_flag == 1) {&exception::error("different line length detected: $id_key", $_[0]);}
			if ($error_flag == 2) {&exception::error("different EOL code detected: $id_key", $_[0]);}
			if ($row_width and $row_width != length($line)) {$error_flag = 1;}
			if ($row_bytes and $row_bytes != $row_bytes1) {$error_flag = 2;}
			
			# 配列長を追加
			$seq_length += length($line);
			if (!$row_width) {$row_width = length($line);}
			if (!$row_bytes) {$row_bytes = $row_bytes1;}
		}
		
		# ID行およびファイル末の処理
		if ($line =~ /^>/ or eof(FASTA)) {
			if ($id_key) {
				# 直前のデータをインデックスハッシュに登録
				$fasta_index{$id_key} = {
					"id_start" => $id_start,
					"seq_length" => $seq_length,
					"seq_start" => $seq_start,
					"row_width" => $row_width,
					"row_bytes" => $row_bytes
				};
				
				# 直前のデータをインデックスファイルに書き込む
				print FASTA_INDEX "$id_key\t$seq_length\t$seq_start\t$row_width\t$row_bytes\n";
			}
			
			# IDの最初の空白文字の直前までをハッシュキーとして登録
			($id_key) = split(/\s/, substr($line, 1));
			
			# 配列開始点を更新
			$seq_start = tell(FASTA);
			
			# ID開始点を更新
			$id_start = $seq_start - $row_bytes1;
			
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
	
	# インデックスハッシュのリファレンスを返す
	return(\%fasta_index);
}

# bedデータ整理 common::modify_bed(bedデータリストリファレンス, 領域番号, 領域名, 出力遺伝子タイプ, 出力バリアント数, 出力データ形式)
sub modify_bed {
	# 変数を宣言
	my %type = ("red" => substr($_[3], 0, 1), "yellow" => substr($_[3], 1, 1), "blue" => substr($_[3], 2, 1), "." => 1);
	my @loci = ();
	
	# データを順に整理
	foreach my $bed (@{$_[0]}) {
		# 未指定の遺伝子タイプのデータを除外
		if (!$type{$bed->[8]}) {next;}
		
		# ブロック情報をリストに変換
		$bed->[10] = [split(/,/, $bed->[10])];
		$bed->[11] = [split(/,/, $bed->[11])];
		
		# ブロック情報を絶対位置に変換して最後列に追加
		my @blocks = ();
		for (my $i = 0;$i < $bed->[9];$i++) {$blocks[$i] = [$bed->[1] + $bed->[11]->[$i], $bed->[1] + $bed->[10]->[$i] + $bed->[11]->[$i]];}
		
		# 変数を宣言
		my $basal_locus = "";
		
		# 既出領域と比較
		for (my $i = 0;$i < @loci;$i++) {
			# 方向が異なる場合を除外
			if ($bed->[5] ne $loci[$i]->[0]->[3]) {next;}
			
			# 領域がオーバーラップしない場合を除外
			if ($bed->[2] <= $loci[$i]->[0]->[1] or $loci[$i]->[0]->[2] <= $bed->[1]) {next;}
			
			# ブロックがオーバーラップしない場合を除外
			if (!List::Util::first {my $block = $_;List::Util::first {$_->[0] < $block->[1] and $block->[0] < $_->[1]} @{$loci[$i]->[0]->[4]}} @blocks) {next;}
			
			# 既存アイソフォームから完全に一致するものを探索
			my $synonym = List::Util::first {join(",", ($bed->[1], @{$bed->[10]}, @{$bed->[11]})) eq join(",", ($_->[1], @{$_->[10]}, @{$_->[11]}))} @{$loci[$i]->[1]};
			
			# 既存アイソフォームから完全に一致する場合を除外
			if ($synonym) {
				if ($bed->[4] > $synonym->[4]) {($synonym->[3], $synonym->[4]) = ($bed->[3], $bed->[4]);}
				$basal_locus = 1;
				last;
			}
			
			# オーバーラップする領域を統合
			if ($basal_locus) {
				if ($loci[$i]->[0]->[1] < $basal_locus->[0]->[1]) {$basal_locus->[0]->[1] = $loci[$i]->[0]->[1];}
				if ($loci[$i]->[0]->[2] > $basal_locus->[0]->[2]) {$basal_locus->[0]->[2] = $loci[$i]->[0]->[2];}
				push(@{$basal_locus->[0]->[4]}, @{$loci[$i]->[0]->[4]});
				push(@{$basal_locus->[1]}, @{$loci[$i]->[1]});
				splice(@loci, $i, 1);
				$i--;
			}
			
			# オーバーラップする領域に追加
			else {
				$basal_locus = $loci[$i];
				if ($bed->[1] < $basal_locus->[0]->[1]) {$basal_locus->[0]->[1] = $bed->[1];}
				if ($bed->[2] > $basal_locus->[0]->[2]) {$basal_locus->[0]->[2] = $bed->[2];}
				push(@{$basal_locus->[0]->[4]}, @blocks);
				push(@{$basal_locus->[1]}, $bed);
			}
		}
		
		# 新規領域を登録
		if (!$basal_locus) {push(@loci, [[$bed->[0], $bed->[1], $bed->[2], $bed->[5], \@blocks], [$bed]]);}
	}
	
	# アイソフォーム部分を抽出
	@loci = map {$_->[1]} @loci;
	
	# 各領域のアイソフォームを並べ替え (スコア順 > 領域長順 > 名前順)
	@loci = map {[(sort {$b->[4] <=> $a->[4] or $b->[2] - $b->[1] <=> $a->[2] - $a->[1] or $a->[3] cmp $b->[3]} @{$_})]} @loci;
	
	# 各領域のアイソフォームを指定した個数だけ選抜
	foreach (@loci) {if ($_[4]) {splice(@{$_}, $_[4]);}}
	
	# 各領域のアイソフォームを位置で並べ替え
	@loci = map {[(sort {$a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] or $a->[3] cmp $b->[3]} @{$_})]} @loci;
	
	# 各領域を位置で並べ替え
	@loci = sort {$a->[0]->[1] <=> $b->[0]->[1] or $a->[0]->[2] <=> $b->[0]->[2] or $a->[0]->[5] cmp $b->[0]->[5]} @loci;
	
	# データを出力
	foreach my $locus (@loci) {
		# 領域番号を加算
		$_[1]++;
		
		# アイソフォームを順に処理
		for (my $i = 0;$i < @{$locus};$i++) {
			# bed12形式で出力
			if ($_[5] eq "bed") {
				$locus->[$i]->[3] = "$_[2]$_[1]." . ($i + 1) . ":" . substr($locus->[$i]->[3], index($locus->[$i]->[3], ":") + 1);
				print join("\t", @{$locus->[$i]}[0..9]), "\t", join(",", @{$locus->[$i]->[10]}), "\t", join(",", @{$locus->[$i]->[11]}), "\n";
			}
			
			# gtf形式で出力
			elsif ($_[5] eq "gtf") {
				# 開始点をgtf形式に修正
				$locus->[$i]->[1]++;
				
				# 遺伝子タイプを修正
				if ($locus->[$i]->[8] eq "red") {$locus->[$i]->[8] = "pseudogene";}
				elsif ($locus->[$i]->[8] eq "yellow") {$locus->[$i]->[8] = "truncated";}
				elsif ($locus->[$i]->[8] eq "blue") {$locus->[$i]->[8] = "protein_coding";}
				elsif ($locus->[$i]->[8] eq ".") {$locus->[$i]->[8] = "uncharacterized";}
				
				# 領域を示す行を出力
				print "$locus->[$i]->[0]\tfate\ttranscript\t$locus->[$i]->[1]\t$locus->[$i]->[2]\t$locus->[$i]->[4]\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $locus->[$i]->[8], '";', "\n";
				
				# 逆鎖の場合は個々のブロックを逆順に並べ替え
				if ($locus->[$i]->[5] eq "-") {$locus->[$i]->[10] = [reverse(@{$locus->[$i]->[10]})];$locus->[$i]->[11] = [reverse(@{$locus->[$i]->[11]})];}
				
				# 変数を宣言
				my $sign = $locus->[$i]->[5] . 1;
				
				# 個々のヒットを示す行を出力
				for (my $j = 0;$j < $locus->[$i]->[9];$j++) {
					# エキソン行を出力
					print "$locus->[$i]->[0]\tfate\texon\t", $locus->[$i]->[1] + $locus->[$i]->[11]->[$j], "\t", $locus->[$i]->[1] + $locus->[$i]->[10]->[$j] + $locus->[$i]->[11]->[$j] - 1, "\t.\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; exon_number "', $j + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $locus->[$i]->[8], '"; exon_id "', "$_[2]$_[1].", $i + 1, ".", $j + 1, '";', "\n";
					
					# コーディング領域以外を除外
					if ($locus->[$i]->[8] ne "protein_coding") {next;}
					
					# コーディング領域行を出力
					print "$locus->[$i]->[0]\tfate\tCDS\t", $locus->[$i]->[1] + $locus->[$i]->[11]->[$j] + 3 * 0 ** ($locus->[$i]->[9] - $j + $sign), "\t", $locus->[$i]->[1] + $locus->[$i]->[10]->[$j] + $locus->[$i]->[11]->[$j] - 1 - 3 * 0 ** ($locus->[$i]->[9] - $j - $sign), "\t.\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; exon_number "', $j + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $locus->[$i]->[8], '"; protein_id "', "$_[2]$_[1].", $i + 1, 'p";', "\n";
				}
			}
		}
	}
	return(1);
}

# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
