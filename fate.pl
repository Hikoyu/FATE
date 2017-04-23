#!/usr/bin/env perl

# Copyright (c) 2016-2017 Hikoyu Suzuki
# This software is released under the MIT License.

use strict;
use warnings;
use Getopt::Std;

# ソフトウェアを定義
### 編集範囲 開始 ###
my $software = "fate.pl";	# ソフトウェアの名前
my $version = "ver.1.3.2";	# ソフトウェアのバージョン
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
	$option{"g STR "} = "Gene prediction engine <exonerate|genewise>";
	$option{"h STR "} = "Homology search engine <tblastn|blastn|dc-megablast|megablast> [blastn]";
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
	if ($opt{"g"} and $opt{"g"} ne "exonerate" and $opt{"g"} ne "genewise") {&exception::error("unknown engine specified: -g $opt{g}");}
	if ($opt{"h"} ne "tblastn" and $opt{"h"} ne "blastn" and $opt{"h"} ne "dc-megablast" and $opt{"h"} ne "megablast") {&exception::error("unknown engine specified: -h $opt{h}");}
	if ($opt{"s"} and !$opt{"g"}) {&exception::caution("-s ignored under -g unspecified");}
	if ($opt{"5"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -5 $opt{5}");}
	if ($opt{"3"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -3 $opt{3}");}
	if ($opt{"i"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -i $opt{i}");}
	if ($opt{"o"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -o $opt{o}");}
	if ($opt{"l"} !~ /^\d+$/) {&exception::error("specify INT >= 0: -l $opt{l}");}
	if ($opt{"c"} !~ /^\d+$|^\d+\.\d+$|^\d+[eE]-?\d+$|^\d+\.\d+[eE]-?\d+$/ or $opt{"c"} > 1) {&exception::error("specify NUM 0-1: -c $opt{c}");}
	
	# 依存するソフトウェアがインストールされているか確認
	if ($opt{"g"} and !`which $opt{"g"}`) {&exception::error("$opt{g} not installed");}
	
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
	my $search_engine = $opt{"h"} =~ /^tblastn/ ? "tblastn" : "blastn -task $opt{h}";
	my $homology_search = "$search_engine -num_threads $opt{p} -outfmt 6 -soft_masking true -db $genome_file";
	my $gene_prediction = $opt{"g"};
	if ($opt{"g"} and $opt{"g"} eq "exonerate") {
		$gene_prediction .= " -m protein2genome -V 0 --showtargetgff T --showalignment F --showvulgar F";
		$gene_prediction .= $opt{"s"} ? " --forcegtag T" : "";
	}
	elsif ($opt{"g"} and $opt{"g"} eq "genewise") {
		$gene_prediction .= " -divide '#' -pseudo -sum -gff -silent -quiet";
		$gene_prediction .= $opt{"s"} ? "" : " -nosplice_gtag";
	}
	$homology_search .= `which tee` ? " | tee fate_search_$opt{h}.out" : "";
	
	# 個々の配列データを保存しておくディレクトリを作成 (-g指定時)
	if ($opt{"g"}) {
		if (!-d "queries") {mkdir("queries") or &exception::error("failed to make directory: queries");}
		if (!-d "loci") {mkdir("loci") or &exception::error("failed to make directory: loci");}
	}
	
	# 変数を宣言
	my @pid = ();
	my @pipe = ();
	
	# プロセス間通信のパイプを作成
	my $order = IO::Pipe->new;
	my $output = IO::Pipe->new;
	
	# 指定したプロセス数で並列処理
	for (my $pnum = 0;$pnum < $opt{"p"};$pnum++) {
		# プロセス間通信のパイプを作成
		push(@pipe, my $input = IO::Pipe->new);
		
		# プロセス分岐
		$pid[$pnum] = fork;
		
		# プロセス分岐に失敗した場合
		if (!defined($pid[$pnum])) {&exception::error("failed to fork process");}
		
		## ここから子プロセスAの処理 ##
		if (!$pid[$pnum]) {
			# パイプを開く
			$order->writer;
			$output->writer;
			$input->reader;
			
			# 子プロセスCにプロセス番号を送信
			syswrite($order, "$pnum\n");
			
			# 遺伝子構造予測を実行 (-g指定時)
			if ($opt{"g"}) {
				# ゲノム配列のfastaファイルを開く
				open(GENOME, "<", $genome_file) or &exception::error("failed to open file: $genome_file");
				
				# 子プロセスCから領域データを受信し、遺伝子構造予測を実行
				while (<$input>) {
					# 改行コードを除去
					chomp;
					
					# タブ文字でデータを分割
					my @col = split(/\t/);
					
					# 1列のデータを受信した場合
					if (@col == 1) {last;}
					
					# マスクブロックサイズとマスクブロック開始点リストを取得
					my @mask_block_size = $col[6] ? split(/,/, $col[7]) : ();
					my @mask_block_start = $col[6] ? split(/,/, $col[8]) : ();
					
					# フランキング配列長に合わせて開始点と終了点を修正
					my $locus_start = $col[5] > 0 ? $col[1] - $opt{"5"} : $col[1] - $opt{"3"};
					my $locus_end = $col[5] > 0 ? $col[2] + $opt{"3"} : $col[2] + $opt{"5"};
					if ($locus_start < 0) {$locus_start = 0;}
					if ($locus_end > $genome_faidx->{$col[0]}->{"seq_length"}) {$locus_end = $genome_faidx->{$col[0]}->{"seq_length"};}
					
					# 取り出す領域の位置情報を算出
					my $start_point = int($locus_start / $genome_faidx->{$col[0]}->{"row_width"}) * $genome_faidx->{$col[0]}->{"row_bytes"} + $locus_start % $genome_faidx->{$col[0]}->{"row_width"};
					my $end_point = int($locus_end / $genome_faidx->{$col[0]}->{"row_width"}) * $genome_faidx->{$col[0]}->{"row_bytes"} + $locus_end % $genome_faidx->{$col[0]}->{"row_width"};
					
					# 領域の塩基配列を取得
					seek(GENOME, $genome_faidx->{$col[0]}->{"seq_start"} + $start_point, 0);
					read(GENOME, my $locus_seq, $end_point - $start_point);
					
					# 改行コードを除去
					$locus_seq =~ s/\n|\r//g;
					
					# 配列を大文字に変換
					$locus_seq = uc($locus_seq);
					
					# 配列をマスク
					for (my $i = 0;$i < $col[6];$i++) {substr($locus_seq, $mask_block_start[$i] - $locus_start, $mask_block_size[$i]) = "N" x $mask_block_size[$i];}
					
					# 相補鎖に変換
					if ($col[5] < 0) {&common::complementary($locus_seq);}
					
					# 領域開始点を表示形式に合わせる
					$locus_start++;
					
					# 領域名を定義
					my $locus_name = "$col[0]:$locus_start-$locus_end(";
					$locus_name .= $col[5] > 0 ? "+)" : "-)";
					
					# 領域の配列ファイルを作成
					open(LOCUS, ">", "loci/locus$pnum.fa") or &exception::error("failed to make file: loci/locus$pnum.fa");
					
					# 取得した配列をファイルに出力
					print LOCUS ">$locus_name\n";
					for (my $pos = 0;$pos < length($locus_seq);$pos += 60) {print LOCUS substr($locus_seq, $pos, 60), "\n";}
					
					# 領域の配列ファイルを閉じる
					close(LOCUS);
					
					# 変数を宣言
					my @genes = ();
					my $query_name = $col[3];
					my $query_start = 0;
					my $query_end = 0;
					my $summary_flag = 0;
					
					# 遺伝子構造予測を実行
					$query_name =~ s/\|/\\\|/g;
					my $args = $opt{"g"} eq "exonerate" ? "-q queries/$query_name.fa -t loci/locus$pnum.fa" : $opt{"g"} eq "genewise" ? "queries/$query_name.fa loci/locus$pnum.fa" : "";
					open(PREDICT, "-|", "$gene_prediction $args 2>/dev/null") or &exception::error("failed to execute gene prediction: $col[3] vs $locus_name");
					
					# 領域開始点をbed形式に戻す
					$locus_start--;
					
					# 遺伝子構造予測結果を1行ずつ読み込んで処理
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
					
					# 遺伝子構造予測を終了
					close(PREDICT);
					
					# 翻訳領域を推定
					foreach my $gene (@genes) {
						# 変数を宣言
						my $upstream_truncation = 0;
						my $downstream_truncation = 0;
						my $completeness = 0;
						
						# Nを含まない上流配列を取得し、その長さが指定値より短い場合は5'側truncateとする
						my $upstream_seq = substr(substr($locus_seq, 0, $gene->[1]), -$opt{"5"});
						my $upstream_offset = $col[5] > 0 ? $col[1] - $locus_start - $gene->[1] : $locus_end - $gene->[1] - $col[2];
						if ($upstream_offset < 0) {$upstream_offset = 0;}
						$upstream_seq =~ s/^.*N//;
						if (length($upstream_seq) < $opt{"5"} - $upstream_offset) {$upstream_truncation = 1;}
						
						# Nを含まない下流配列を取得し、その長さが指定値より短い場合は3'側truncateとする
						my $downstream_seq = substr($locus_seq, $gene->[2], $opt{"3"});
						my $downstream_offset = $col[5] > 0 ? $locus_start + $gene->[2] - $col[2] : $col[1] - $locus_end + $gene->[2];
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
						
						# 領域開始点・終了点を修正
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
						($gene->[1], $gene->[2]) = ($col[5] > 0 ? $locus_start + $gene->[1] : $locus_end - $gene->[2], $col[5] > 0 ? $locus_start + $gene->[2] : $locus_end - $gene->[1]);
						($gene->[6], $gene->[7]) = ($gene->[1], $gene->[2]);
						$gene->[5] = $col[5];
						if ($col[5] < 0) {
							$gene->[10] = [reverse(@{$gene->[10]})];
							$gene->[11] = [reverse(@{$gene->[11]})];
							for (my $i = 0;$i < @{$gene->[11]};$i++) {$gene->[11]->[$i] = $gene->[2] - $gene->[1] - $gene->[10]->[$i] - $gene->[11]->[$i];}
						}
						
						# 子プロセスBにbed12形式でデータを送信
						syswrite($output, join("\t", (@{$gene}[0..9], join(",", @{$gene->[10]}), join(",", @{$gene->[11]}))) . "\n");
					}
					
					# 子プロセスBにデータ区切りを送信
					syswrite($output, "$locus_name\n");
					
					# 子プロセスCにプロセス番号を送信
					syswrite($order, "$pnum\n");
				}
				
				# ゲノム配列のfastaファイルを閉じる
				close(GENOME);
				
				# 子プロセスBにデータ転送完了を送信
				syswrite($output, "//\n");
			}
			
			# 変数を宣言
			my @loci = ();
			
			# 子プロセスBからbedデータを受信しながらアイソフォームを決定
			while (<$input>) {
				# 改行コードを除去
				chomp;
				
				# タブ文字でデータを分割
				my @col = split(/\t/);
				
				# 1列のデータを受信した場合
				if (@col == 1) {
					# データ転送完了を受信した場合
					if ($col[0] eq "//") {last;}
					
					# 遺伝子を分離
					my $genes = &common::define_isoforms(\@loci, $opt{"v"}, $opt{"t"});
					
					# 親プロセスにデータを送信
					foreach my $gene (@{$genes}) {
						syswrite($output, "$col[0]\t" . join("\t", @{$gene->[0]}) . "\n");
						foreach my $isoform (@{$gene->[1]}) {syswrite($output, "$col[0]\t" . join("\t", @{$isoform}[0..9]) . "\t" . join(",", @{$isoform->[10]}) . "\t" . join(",", @{$isoform->[11]}) . "\n");}
					}
					
					# 親プロセスに通し番号を送信
					syswrite($output, "$col[0]\n");
					
					# リストをリセット
					@loci = ();
					
					# 子プロセスBにプロセス番号を送信
					syswrite($order, "$pnum\n");
					next;
				}
				
				# ブロック情報をリストに変換
				$col[10] = [split(/,/, $col[10])];
				$col[11] = [split(/,/, $col[11])];
				
				# データをリストに追加
				push(@loci, \@col);
			}
			
			# プロセスを終了
			exit 0;
		}
		## ここまで子プロセスAの処理 ##
	}
	
	# プロセス間通信のパイプを作成
	my $signal = IO::Pipe->new;
	
	# プロセス分岐
	$pid[$opt{"p"}] = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid[$opt{"p"}])) {&exception::error("failed to fork process");}
	
	## ここから子プロセスBの処理 ##
	if (!$pid[$opt{"p"}]) {
		# パイプを開く
		$output->reader;
		
		# 変数を宣言
		my %loci = ();
		my $fin_loci = 0;
		my $fin_procs = 0;
		
		# 子プロセスCからアセンブルID数を受信
		my $num_loci = <$output>;
		
		# 子プロセスAまたはCからデータを受信してハッシュに登録
		print STDERR $opt{"g"} ? "Running gene prediction...0%" : "Loading assembled loci...";
		while (<$output>) {
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# 1列のデータを受信した場合
			if (@col == 1) {
				# データ転送完了を受信した場合
				if ($col[0] eq "//\n") {
					$fin_procs++;
					if ($opt{"g"} and $fin_procs < $opt{"p"}) {next;}
					else {last;}
				}
				
				# データ区切りを受信した場合
				else {
					$fin_loci++;
					print STDERR "\rRunning gene prediction...", int($fin_loci / $num_loci * 100), "%";
					next;
				}
			}
			
			# bed形式のデータをハッシュに登録
			push(@{$loci{$col[0]}}, $_);
		}
		print STDERR $opt{"g"} ? "\rRunning gene prediction...completed\n" : "completed\n";
		
		# パイプを開く
		$signal->writer;
		
		# 子プロセスCにシグナルを受信
		syswrite($signal, "//\n");
		
		# 変数を宣言
		my $active_pnum = 0;
		my $fin_chr = 0;
		my $num_chr = keys(%loci);
		
		# パイプを開く
		$order->reader;
		foreach (@pipe) {$_->writer;}
		
		# ゲノムデータのID順に処理
		print STDERR "Running isoform definition...0%";
		foreach my $subject (sort {$genome_faidx->{$a}->{"id_order"} <=> $genome_faidx->{$b}->{"id_order"}} keys(%loci)) {
			# 子プロセスAからプロセス番号を受信
			$active_pnum = <$order>;
			
			# アクティブな子プロセスAにbed12形式でデータを送信
			foreach (@{$loci{$subject}}) {syswrite($pipe[$active_pnum], $_);}
			
			# アクティブな子プロセスAに通し番号を送信
			syswrite($pipe[$active_pnum], "$fin_chr\n");
			
			# 通し番号を加算
			$fin_chr++;
			
			# 経過を表示
			print STDERR "\rRunning isoform definition...", int($fin_chr / $num_chr * 100), "%";
		}
		print STDERR "\rRunning isoform definition...completed\n";
		
		# 全ての子プロセスAにデータ転送完了を送信
		foreach (@pipe) {syswrite($_, "//\n");}
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスBの処理 ##
	
	# 相同性検索の出力ファイルを開く (-b指定時)
	if ($opt{"b"}) {open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}");}
	
	# 相同性検索を実行 (-b未指定時)
	else {main::open2(*SEARCH_OUT, *SEARCH_IN, $homology_search) or &exception::error("failed to execute homology search");}
	
	# プロセス間通信のパイプを作成
	my $report = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid)) {&exception::error("failed to fork process");}
	
	## ここから子プロセスCの処理 ##
	if (!$pid) {
		# 相同性検索への入力を閉じる (-b未指定時)
		if (!$opt{"b"}) {close(SEARCH_IN);}
		
		# 変数を宣言
		my %assembled_hits = ();
		my %query_len = ();
		my @blast_hits = ();
		my $last_query = "";
		my $last_subject = "";
		my $num_assembly = 0;
		
		# パイプを開く
		$report->reader;
		
		# 相同性検索の出力を読み込みながら処理
		print STDERR $opt{"b"} ? "Loading homology search results..." : "Running homology search...";
		while (<SEARCH_OUT>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# クエリー名またはサブジェクト名が変わった場合
			if ($col[0] ne $last_query or $col[1] ne $last_subject ) {
				# ヒットをアセンブル
				if ($last_query and $last_subject) {
					$assembled_hits{$last_query}{$last_subject} = &assemble(\@blast_hits, $opt{"i"}, $opt{"o"});
					$num_assembly += scalar(map {keys(%{$_->{"assemble"}})} @{$assembled_hits{$last_query}{$last_subject}});
				}
				
				# クエリー名が変わった場合
				if ($col[0] ne $last_query) {
					# クエリー名を更新
					$last_query = $col[0];
					
					# 親プロセスからクエリー配列長を受信して登録
					$query_len{$col[0]} = <$report>;
					
					# 改行コードを除去
					chomp($query_len{$col[0]});
				}
				
				# サブジェクト名を更新
				$last_subject = $col[1];
				
				# リストをリセット
				@blast_hits = ();
			}
			
			# データをリストに登録
			push(@blast_hits, $col[8] < $col[9] ? {"query_start" => $col[6] - 1, "query_end" => $col[7], "locus_start" => $col[8] - 1, "locus_end" => $col[9], "strand" => 1} : {"query_start" => $col[7], "query_end" => $col[6] - 1, "locus_start" => $col[9] - 1, "locus_end" => $col[8], "strand" => -1});
			$blast_hits[-1]->{"assemble"} = {};
			$blast_hits[-1]->{"score"} = $col[11];
			$blast_hits[-1]->{"num_connection"} = 0;
		}
		print STDERR "completed\n";
		
		# 相同性検索の出力を閉じる
		close(SEARCH_OUT);
		
		# 残りのヒットをアセンブル
		if ($last_query and $last_subject) {
			$assembled_hits{$last_query}{$last_subject} = &assemble(\@blast_hits, $opt{"i"}, $opt{"o"});
			$num_assembly += scalar(map {keys(%{$_->{"assemble"}})} @{$assembled_hits{$last_query}{$last_subject}});
		}
		
		# 相同性検索でヒットが得られなかった場合
		else {&exception::caution("no hits found from $opt{h} search");}
		
		# パイプを開く
		$output->writer;
		
		# 子プロセスBにアセンブルID数を送信
		syswrite($output, "$num_assembly\n");
		
		# アクティブな子プロセスAにbed6形式+マスクブロック情報でデータを送信 (-g指定時)
		if ($opt{"g"}) {
			# パイプを開く
			$order->reader;
			foreach (@pipe) {$_->writer;}
			
			# アセンブルデータを順に送信
			foreach my $query (keys(%assembled_hits)) {
				foreach my $subject (keys(%{$assembled_hits{$query}})) {
					foreach my $locus (@{$assembled_hits{$query}{$subject}}) {
						foreach my $assemble (values(%{$locus->{"assemble"}})) {
							syswrite($pipe[<$order>], join("\t", ($subject, $locus->{"locus_start"}, $assemble->{"locus_destination"}, $query, $opt{"h"} =~ /^tblastn/ ? $query_len{$query} * 3 : $query_len{$query}, $locus->{"strand"}, scalar(@{$assemble->{"mask_block_size"}}), join(",", @{$assemble->{"mask_block_size"}}), join(",", @{$assemble->{"mask_block_start"}}))) . "\n");
						}
					}
				}
			}
			
			# 全ての子プロセスAにデータ転送完了を送信
			foreach (@pipe) {syswrite($_, "//\n");}
		}
		
		# 子プロセスBにbed12形式でデータを送信 (-g未指定時)
		else {
			# アセンブルデータを順に送信
			foreach my $query (keys(%assembled_hits)) {
				foreach my $subject (keys(%{$assembled_hits{$query}})) {
					foreach my $locus (@{$assembled_hits{$query}{$subject}}) {
						foreach my $assemble (values(%{$locus->{"assemble"}})) {
							syswrite($output, join("\t", ($subject, $locus->{"locus_start"}, $assemble->{"locus_destination"}, $query, $assemble->{"total_score"}, $locus->{"strand"}, $locus->{"locus_start"}, $assemble->{"locus_destination"}, ".", scalar(@{$assemble->{"block_size"}}), join(",", @{$assemble->{"block_size"}}), join(",", map {$_ - $locus->{"locus_start"}} @{$assemble->{"block_start"}}))) . "\n");
						}
					}
				}
			}
			
			# 子プロセスBにデータ転送完了を送信
			syswrite($output, "//\n");
		}
		
		# パイプを開く
		$signal->reader;
		
		# 子プロセスBからシグナルを受信
		<$signal>;
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスCの処理 ##
	
	# 相同性検索の出力を閉じる
	close(SEARCH_OUT);
	
	# 変数を宣言
	my $query_title = "";
	my $query_seq = "";
	
	# パイプを開く
	$report->writer;
	
	# クエリー配列を読み込みながら処理
	while (my $line = <>) {
		# 改行コードを除去
		$line =~ s/\n|\r//g;
		
		# 配列行の処理
		$query_seq .= $line =~ /^>/ ? "" : uc($line);
		
		# ID行の処理
		if ($line =~ /^>/ or eof) {
			# 配列データを処理
			if ($query_title) {
				# 子プロセスCにクエリー配列長を送信
				syswrite($report, length($query_seq) . "\n");
				
				# クエリー配列をファイルに出力 (-g指定時)
				if ($opt{"g"}) {
					# 変数を宣言
					my $output_seq = $opt{"h"} =~ /^tblastn/ ? $query_seq : &common::translate($query_seq, 0);
					
					# 配列からアスタリスクを除去
					$output_seq =~ s/\*//g;
					
					# 出力ファイルを作成
					open(QUERY, ">", "queries/$query_title.fa") or &exception::error("failed to make file: queries/$query_title.fa");
					
					# クエリー配列をfasta形式でファイルに出力
					print QUERY ">$query_title\n";
					for (my $pos = 0;$pos < length($output_seq);$pos += 60) {print QUERY substr($output_seq, $pos, 60), "\n";}
					
					# 出力ファイルを閉じる
					close(QUERY);
				}
				
				# クエリー配列を相同性検索に入力 (-b未指定時)
				if (!$opt{"b"}) {
					print SEARCH_IN ">$query_title\n";
					for (my $pos = 0;$pos < length($query_seq);$pos += 60) {print SEARCH_IN substr($query_seq, $pos, 60), "\n";}
				}
			}
			
			# ID行の最初の空白文字の前までを配列名として登録
			($query_title) = split(/\s/, substr($line, 1));
			
			# 配列をリセット
			$query_seq = "";
		}
	}
	
	# 相同性検索への入力を閉じる (-b未指定時)
	if (!$opt{"b"}) {close(SEARCH_IN);}
	
	# 子プロセスCを刈り取る
	waitpid($pid, 0);
	if ($?) {&exception::error("process abnormally exited");}
	
	# 変数を宣言
	my @loci = ();
	my @fin_loci = ();
	my $fin_chr = 0;
	my $gene_num = 0;
	
	# パイプを開く
	$output->reader;
	
	# データを受信しながら処理
	while (<$output>) {
		# 改行コードを除去
		chomp;
		
		# タブ文字でデータを分割
		my @col = split(/\t/);
		
		# bed12形式のデータを受信した場合
		if (@col > 12) {
			# ブロック情報をリストに変換
			$col[11] = [split(/,/, $col[11])];
			$col[12] = [split(/,/, $col[12])];
			
			# アイソフォームとしてデータを追加
			push(@{$loci[$col[0] - $fin_chr]->[-1]}, [@col[1..12]]);
		}
		
		# bed4形式のデータを受信した場合
		elsif (@col > 4) {
			# 遺伝子としてデータを追加
			push(@{$loci[$col[0] - $fin_chr]}, [[@col[1..4]]]);
		}
		
		# 通し番号を受信した場合
		else {
			# 受信した通し番号の領域のデータ受信が完了したことを登録
			$fin_loci[$col[0] - $fin_chr] = 1;
			
			# 通し番号順に領域のデータ受信が完了していればデータを出力
			while ($fin_loci[0]) {
				&common::output($loci[0], $gene_num, $opt{"n"}, $opt{"f"});
				shift(@loci);
				shift(@fin_loci);
				$fin_chr++;
			}
		}
	}
	
	# 子プロセスA、Bを刈り取る
	foreach (@pid) {
		waitpid($_, 0);
		if ($?) {&exception::error("process abnormally exited");}
	}
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
			if ($sorted_hits[$k]->{"locus_start"} - $sorted_hits[$i]->{"locus_end"} > $_[1]) {last;}
			
			# 方向が異なるヒットを除外
			if (!($sorted_hits[$k]->{"strand"} + $sorted_hits[$i]->{"strand"})) {next;}
			
			# 指定値を超えるクエリーオーバーラップのヒットを除外
			if (abs($sorted_hits[$k]->{"query_start"} - $sorted_hits[$i]->{"query_end"}) > $_[2]) {next;}
			
			# クエリー開始点が前進しないヒットを除外
			if ($sorted_hits[$k]->{"query_start"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_start"} * $sorted_hits[$i]->{"strand"}) {next;}
			
			# クエリー終了点が前進しないヒットを除外
			if ($sorted_hits[$k]->{"query_end"} * $sorted_hits[$k]->{"strand"} <= $sorted_hits[$i]->{"query_end"} * $sorted_hits[$i]->{"strand"}) {next;}
			
			# マスクするヒットを決定
			my @mask_block = grep {$_->{"strand"} + $sorted_hits[$i]->{"strand"}} @sorted_hits[$i + 1..$k - 1];
			my @mask_block_size = map {$_->{"locus_end"} - $_->{"locus_start"}} @mask_block;
			my @mask_block_start = map {$_->{"locus_start"}} @mask_block;
			
			# 各アセンブルIDについて処理
			foreach (keys(%{$sorted_hits[$k]->{"assemble"}})) {
				# 総スコアが更新されない同一IDのアセンブルを除外
				if (exists($sorted_hits[$i]->{"assemble"}->{$_}) and $sorted_hits[$i]->{"assemble"}->{$_}->{"total_score"} >= $sorted_hits[$i]->{"score"} + $sorted_hits[$k]->{"assemble"}->{$_}->{"total_score"}) {next;}
				
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
		
		# 末尾ヒットの場合
		if (!%{$sorted_hits[$i]->{"assemble"}}) {
			# データを登録
			$sorted_hits[$i]->{"assemble"}->{$assemble_id} = {
				"total_score" => $sorted_hits[$i]->{"score"},
				"locus_destination" => $sorted_hits[$i]->{"locus_end"},
				"block_size" => [$sorted_hits[$i]->{"locus_end"} - $sorted_hits[$i]->{"locus_start"}],
				"block_start" => [$sorted_hits[$i]->{"locus_start"}],
				"mask_block_size" => [],
				"mask_block_start" => []
			};
			
			# アセンブルIDを更新
			$assemble_id++;
		}
	}
	
	# アセンブルデータリストのリファレンスを返す
	return([grep {!$_->{"num_connection"}} @sorted_hits]);
}

## ここからfilterコマンドのパッケージ ##
package filter;

# コマンドとオプションを定義
sub define {
	$note = "Filter already annotated protein-coding genes under specified conditions.";
	$usage = "<genome.fa> <STDIN | in1.bed> [in2.bed ...] [> out.bed | > out.gtf]";
	$option{"d PATH "} = "Path to gene or protein database file (fasta format)";
	$option{"h STR "} = "Homology search engine <blastx|blastn|dc-megablast|megablast> [blastx]";
	$option{"k STR "} = "Keywords for filtering (AND[&], OR[;], BUT[!])";
	$option{"r INT "} = "Cutoff rank of hits <1->";
	return(1);
}

# コマンド本体
sub body {
	# 指定されたオプションを確認
	if ($opt{"b"} and $opt{"d"}) {&exception::error("options incompatible: -b and -d");}
	if ($opt{"h"} ne "blastx" and $opt{"h"} ne "blastn" and $opt{"h"} ne "dc-megablast" and $opt{"h"} ne "megablast") {&exception::error("unknown engine specified: -h $opt{h}");}
	if (defined($opt{"r"}) and ($opt{"r"} !~ /^\d+$/ or $opt{"r"} < 1)) {&exception::error("specify INT >= 1: -r $opt{r}");}
	
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
	
	# 参照データファイルを確認 (-d指定時)
	if ($opt{"d"}) {
		# 参照配列のfastaファイルを検索
		&common::find_db($opt{"d"});
		
		# 参照配列のblastデータベースを確認
		if ($opt{"h"} =~ /^blastx/) {&common::check_blastdb($opt{"d"}, 1);} else {&common::check_blastdb($opt{"d"});}
	}
	
	# 処理を定義
	my $search_engine = $opt{"h"} =~ /^blastx/ ? "blastx" : "blastn -task $opt{h}";
	my $homology_search = "$search_engine -num_threads $opt{p} -outfmt '6 std salltitles' -soft_masking true -strand plus";
	$homology_search .= $opt{"d"} ? " -db $opt{d}" : "";
	$homology_search .= `which tee` ? " | tee fate_filter_$opt{h}.out" : "";

	# 変数を宣言
	my @pid = ();
	my @pipe = ();
	
	# プロセス間通信のパイプを作成
	my $order = IO::Pipe->new;
	my $output = IO::Pipe->new;
	
	# 指定したプロセス数で並列処理
	for (my $pnum = 0;$pnum < $opt{"p"};$pnum++) {
		# プロセス間通信のパイプを作成
		push(@pipe, my $input = IO::Pipe->new);
		
		# プロセス分岐
		$pid[$pnum] = fork;
		
		# プロセス分岐に失敗した場合
		if (!defined($pid[$pnum])) {&exception::error("failed to fork process");}
		
		## ここから子プロセスAの処理 ##
		if (!$pid[$pnum]) {
			# 変数を宣言
			my %blast_hits = ();
			
			# パイプを開く
			$order->writer;
			$output->writer;
			$input->reader;
			
			# 子プロセスCにプロセス番号を送信
			syswrite($order, "$pnum\n");
			
			# 相同性検索の結果を集計 (-bまたは-d指定時)
			if ($opt{"b"} or $opt{"d"}) {
				# 子プロセスCからデータを受信しながら処理
				while (<$input>) {
					# 改行コードを除去
					chomp;
					
					# タブ文字でデータを分割
					my @col = split(/\t/);
					
					# 1列のデータを受信した場合
					if (@col == 1) {
						# データ転送完了を受信した場合
						if ($col[0] eq "//") {last;}
						
						# データ区切りを受信した場合
						else {
							# 変数を宣言
							my %score = ();
							
							# 各ヒットについて処理
							foreach (keys(%blast_hits)) {
								# 変数を宣言
								my @query_array = ();
								my @target_array = ();
								my $total_score = 0;
								my $total_length = 0;
								
								# スコアを算出
								foreach (@{$blast_hits{$_}}) {
									$total_score += $_->{"score"};
									$total_length += $_->{"query_end"} - $_->{"query_start"} + $_->{"locus_end"} - $_->{"locus_start"};
									@query_array[$_->{"query_start"}..$_->{"query_end"} - 1] = (1) x ($_->{"query_end"} - $_->{"query_start"});
									@target_array[$_->{"locus_start"}..$_->{"locus_end"} - 1] = (1) x ($_->{"locus_end"} - $_->{"locus_start"});
								}
								$score{$_} = $total_score / $total_length * scalar(grep {$_} (@query_array, @target_array));
							}
							
							# ヒットをスコア順に並べ替え
							my @subjects = sort {$score{$b} <=> $score{$a}} keys(%score);
							
							# 指定値より下位のヒットを削除
							if ($opt{"r"}) {splice(@subjects, $opt{"r"});}
							
							# 条件を満たした場合は子プロセスBに領域名と真を送信
							if (grep {&common::keyword_search($_, $opt{"k"})} @subjects) {syswrite($output, "$col[0]\t1\n");}
							
							# 条件を満たさなかった場合は子プロセスBに領域名と偽を送信
							else {syswrite($output, "$col[0]\t0\n");}
							
							# ハッシュをリセット
							%blast_hits = ();
							
							# 子プロセスCにプロセス番号を送信
							syswrite($order, "$pnum\n");
							next;
						}
					}
					
					# ヒット名を編集
					if ($col[12]) {$col[1] .= " $col[12]";}
					
					# ヒットをハッシュに登録
					push(@{$blast_hits{$col[1]}}, {"query_start" => $col[6] - 1, "query_end" => $col[7], "locus_start" => $col[8] - 1, "locus_end" => $col[9], "score" => $col[11]});
				}
				
				# 子プロセスBにデータ転送完了を送信
				syswrite($output, "//\n");
			}
			
			# 変数を宣言
			my @loci = ();
			
			# 子プロセスBからbedデータを受信しながらアイソフォームを決定
			while (<$input>) {
				# 改行コードを除去
				chomp;
				
				# タブ文字でデータを分割
				my @col = split(/\t/);
				
				# 1列のデータを受信した場合
				if (@col == 1) {
					# データ転送完了を受信した場合
					if ($col[0] eq "//") {last;}
					
					# 遺伝子を分離
					my $genes = &common::define_isoforms(\@loci, $opt{"v"}, $opt{"t"});
					
					# 親プロセスにデータを送信
					foreach my $gene (@{$genes}) {
						syswrite($output, "$col[0]\t" . join("\t", @{$gene->[0]}) . "\n");
						foreach my $isoform (@{$gene->[1]}) {syswrite($output, "$col[0]\t" . join("\t", @{$isoform}[0..9]) . "\t" . join(",", @{$isoform->[10]}) . "\t" . join(",", @{$isoform->[11]}) . "\n");}
					}
					
					# 親プロセスに通し番号を送信
					syswrite($output, "$col[0]\n");
					
					# リストをリセット
					@loci = ();
					
					# 子プロセスBにプロセス番号を送信
					syswrite($order, "$pnum\n");
					next;
				}
				
				# ブロック情報をリストに変換
				$col[10] = [split(/,/, $col[10])];
				$col[11] = [split(/,/, $col[11])];
				
				# データをリストに追加
				push(@loci, \@col);
			}
			
			# プロセスを終了
			exit 0;
		}
		## ここまで子プロセスAの処理 ##
	}
	
	# プロセス間通信のパイプを作成
	my $signal = IO::Pipe->new;
	
	# プロセス分岐
	$pid[$opt{"p"}] = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid[$opt{"p"}])) {&exception::error("failed to fork process");}
	
	## ここから子プロセスBの処理 ##
	if (!$pid[$opt{"p"}]) {
		# 変数を宣言
		my %bed = ();
		my %loci = ();
		my $fin_procs = 0;
		
		# パイプを開く
		$output->reader;
		
		# 子プロセスAとCからのデータを受信してハッシュに登録
		while (<$output>) {
			# 改行コードを除去
			chomp;
			
			# タブ文字でデータを分割
			my @col = split(/\t/);
			
			# 1列のデータを受信した場合
			if (@col == 1) {
				$fin_procs++;
				if (($opt{"b"} or $opt{"d"}) and $fin_procs < $opt{"p"}) {next;}
				else {last;}
			}
			
			# 2列のデータを受信した場合
			if (@col == 2) {
				# 真を受信した場合はデータを登録
				if ($col[1]) {push(@{$loci{$bed{$col[0]}->[0]}}, $bed{$col[0]});}
				
				# データを削除
				delete($bed{$col[0]});
				next;
			}
			
			# bed形式のデータを一時保存
			$bed{$col[3]} = \@col;
		}
		
		# パイプを開く
		$signal->writer;
		
		# 子プロセスCにシグナルを受信
		syswrite($signal, "//\n");
		
		# 変数を宣言
		my $active_pnum = 0;
		my $fin_chr = 0;
		my $num_chr = keys(%loci);
		
		# パイプを開く
		$order->reader;
		foreach (@pipe) {$_->writer;}
		
		# ゲノムデータのID順に処理
		print STDERR "Running isoform definition...0%";
		foreach my $subject (sort {$genome_faidx->{$a}->{"id_order"} <=> $genome_faidx->{$b}->{"id_order"}} keys(%loci)) {
			# 子プロセスAからプロセス番号を受信
			$active_pnum = <$order>;
			
			# アクティブな子プロセスAにbed12形式でデータを送信
			foreach (@{$loci{$subject}}) {syswrite($pipe[$active_pnum], join("\t", @{$_}) . "\n");}
			
			# アクティブな子プロセスAに通し番号を送信
			syswrite($pipe[$active_pnum], "$fin_chr\n");
			
			# 通し番号を加算
			$fin_chr++;
			
			# 経過を表示
			print STDERR "\rRunning isoform definition...", int($fin_chr / $num_chr * 100), "%";
		}
		print STDERR "\rRunning isoform definition...completed\n";
		
		# 全ての子プロセスAにデータ転送完了を送信
		foreach (@pipe) {syswrite($_, "//\n");}
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスBの処理 ##
	
	# 相同性検索の出力ファイルを開く (-b指定時)
	if ($opt{"b"}) {open(SEARCH_OUT, "<", $opt{"b"}) or &exception::error("failed to open file: $opt{b}");}
	
	# 相同性検索を実行 (-d指定時)
	if ($opt{"d"}) {main::open2(*SEARCH_OUT, *SEARCH_IN, $homology_search) or &exception::error("failed to execute homology search");}
	
	# プロセス間通信のパイプを作成
	my $report = IO::Pipe->new;
	
	# プロセス分岐
	my $pid = fork;
	
	# プロセス分岐に失敗した場合
	if (!defined($pid)) {&exception::error("failed to fork process");}
	
	## ここから子プロセスCの処理 ##
	if (!$pid) {
		# 相同性検索への入力を閉じる (-d指定時)
		if ($opt{"d"}) {close(SEARCH_IN);}
		
		# 変数を宣言
		my @bed = ();
		my $last_query = "";
		my $active_pnum = 0;
		
		# パイプを開く
		$report->reader;
		$output->writer;
		
		# 相同性検索を利用する場合 (-bまたは-d指定時)
		if ($opt{"b"} or $opt{"d"}) {
			# パイプを開く
			$order->reader;
			foreach (@pipe) {$_->writer;}
			
			# 相同性検索の出力を読み込みながら処理
			print STDERR $opt{"b"} ? "Loading homology search results..." : "Running homology search...";
			while (<SEARCH_OUT>) {
				# タブ文字でデータを分割
				my @col = split(/\t/);
				
				# クエリー名を編集
				($col[0]) = split(/::/, $col[0]);
				
				# クエリー名が変わった場合
				if ($col[0] ne $last_query) {
					# 子プロセスAにデータ区切りを送信
					if ($last_query) {syswrite($pipe[$active_pnum], "$last_query\n");}
					
					# クエリー名が一致するまで親プロセスからbedデータを受信
					do {@bed = split(/\t/, <$report>);} until ($bed[3] eq $col[0]);
					
					# 子プロセスBにbedデータを送信
					syswrite($output, join("\t", @bed));
					
					# クエリー名を更新
					$last_query = $col[0];
					
					# 子プロセスAからプロセス番号を受信
					$active_pnum = <$order>;
				}
				
				# アクティブな子プロセスAにデータを送信
				syswrite($pipe[$active_pnum], join("\t", @col));
			}
			print STDERR "completed\n";
			
			# 相同性検索の出力を閉じる
			close(SEARCH_OUT);
			
			# アクティブな子プロセスAに残りのデータとデータ区切りを送信
			if ($last_query) {syswrite($pipe[$active_pnum], "$last_query\n");}
			
			# 相同性検索でヒットが得られなかった場合
			else {&exception::caution("no hits found from $opt{h} search");}
			
			# 全ての子プロセスAにデータ転送完了を送信
			foreach (@pipe) {syswrite($_, "//\n");}
		}
		
		# 相同性検索を利用しない場合 (-b、d未指定時)
		else {
			while (<$report>) {
				# タブ文字でデータを分割
				my @col = split(/\t/);
				
				# 子プロセスBにbedデータと真を送信
				syswrite($output, "$_$col[3]\t1\n");
			}
			
			# 子プロセスBにデータ転送完了を送信
			syswrite($output, "//\n");
		}
		
		# パイプを開く
		$signal->reader;
		
		# 子プロセスBからシグナルを受信
		<$signal>;
		
		# プロセスを終了
		exit 0;
	}
	## ここまで子プロセスCの処理 ##
	
	# 相同性検索の出力を閉じる
	close(SEARCH_OUT);
	
	# パイプを開く
	$report->writer;
	
	# ゲノム配列のfastaファイルを開く
	open(GENOME, "<", $genome_file) or &exception::error("failed to open file: $genome_file");
	
	# bedデータを読み込みながら処理
	while (my $line = <>) {
		# 改行コードを除去
		$line =~ s/\n|\r//g;
		
		# タブ文字でデータを分割
		my @col = split(/\t/, $line);
		
		# 3列未満のデータの場合
		if (@col < 3) {&exception::error("bad format input");}
		
		# 4列目が未定義の場合
		if (!defined($col[3])) {$col[3] = "$opt{n}$.";}
		
		# 5列目が未定義の場合
		if (!defined($col[4])) {$col[4] = 1000;}
		
		# 6列目が未定義の場合
		if (!defined($col[5])) {$col[5] = "+";}
		
		# 7列目が未定義の場合
		if (!defined($col[6])) {$col[6] = $col[1];}
		
		# 8列目が未定義の場合
		if (!defined($col[7])) {$col[7] = $col[2];}
		
		# 9列目が未定義の場合
		if (!defined($col[8])) {$col[8] = ".";}
		
		# 10列目が未定義の場合
		if (!defined($col[9])) {$col[9] = 1;}
		
		# 11列目が未定義の場合
		if (!defined($col[10])) {$col[10] = $col[2] - $col[1];}
		
		# 12列目が未定義の場合
		if (!defined($col[11])) {$col[11] = 0;}
		
		# 向きを数字に変換
		if ($col[5]) {$col[5] = $col[5] . "1";}
		
		# 子プロセスBまたはCにbedデータを送信
		syswrite($report, join("\t", @col) . "\n");
		
		# 相同性検索を実行しない場合 (-d未指定時)
		if (!$opt{"d"}) {next;}
		
		# 取り出す領域の位置情報を算出
		my $start_point = int($col[1] / $genome_faidx->{$col[0]}->{"row_width"}) * $genome_faidx->{$col[0]}->{"row_bytes"} + $col[1] % $genome_faidx->{$col[0]}->{"row_width"};
		my $end_point = int($col[2] / $genome_faidx->{$col[0]}->{"row_width"}) * $genome_faidx->{$col[0]}->{"row_bytes"} + $col[2] % $genome_faidx->{$col[0]}->{"row_width"};
		
		# 領域の塩基配列を取得
		seek(GENOME, $genome_faidx->{$col[0]}->{"seq_start"} + $start_point, 0);
		read(GENOME, my $locus_seq, $end_point - $start_point);
		
		# 改行コードを除去
		$locus_seq =~ s/\n|\r//g;
		
		# 変数を宣言
		my $query_seq = $locus_seq;
		
		# ブロックの領域だけを連結
		my @block_length = split(/,/, $col[10]);
		my @block_pos = split(/,/, $col[11]);
		if ($col[9] == @block_length and $col[9] == @block_pos) {
			$query_seq = "";
			for (my $i = 0;$i < $col[9];$i++) {$query_seq .= substr($locus_seq, $block_pos[$i], $block_length[$i]);}
		}
		
		# 相補鎖に変換
		if ($col[5] and $col[5] < 0) {&common::complementary($query_seq);}
		
		# 配列を相同性検索に入力
		print SEARCH_IN ">$col[3]\n";
		for (my $pos = 0;$pos < length($query_seq);$pos += 60) {print SEARCH_IN substr($query_seq, $pos, 60), "\n";}
	}
	
	# ゲノム配列のfastaファイルを閉じる
	close(GENOME);
	
	# 相同性検索への入力を閉じる (-d指定時)
	if ($opt{"d"}) {close(SEARCH_IN);}
	
	# パイプを閉じる
	$report = undef;
	
	# 子プロセスCを刈り取る
	waitpid($pid, 0);
	if ($?) {&exception::error("process abnormally exited");}
	
	# 変数を宣言
	my @loci = ();
	my @fin_loci = ();
	my $fin_chr = 0;
	my $gene_num = 0;
	
	# パイプを開く
	$output->reader;
	
	# データを受信しながら処理
	while (<$output>) {
		# 改行コードを除去
		chomp;
		
		# タブ文字でデータを分割
		my @col = split(/\t/);
		
		# bed12形式のデータを受信した場合
		if (@col > 12) {
			# ブロック情報をリストに変換
			$col[11] = [split(/,/, $col[11])];
			$col[12] = [split(/,/, $col[12])];
			
			# アイソフォームとしてデータを追加
			push(@{$loci[$col[0] - $fin_chr]->[-1]}, [@col[1..12]]);
		}
		
		# bed4形式のデータを受信した場合
		elsif (@col > 4) {
			# 遺伝子としてデータを追加
			push(@{$loci[$col[0] - $fin_chr]}, [[@col[1..4]]]);
		}
		
		# 通し番号を受信した場合
		else {
			# 受信した通し番号の領域のデータ受信が完了したことを登録
			$fin_loci[$col[0] - $fin_chr] = 1;
			
			# 通し番号順に領域のデータ受信が完了していればデータを出力
			while ($fin_loci[0]) {
				&common::output($loci[0], $gene_num, $opt{"n"}, $opt{"f"});
				shift(@loci);
				shift(@fin_loci);
				$fin_chr++;
			}
		}
	}
	
	# 子プロセスA、Bを刈り取る
	foreach (@pid) {
		waitpid($_, 0);
		if ($?) {&exception::error("process abnormally exited");}
	}
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
		$aa .= exists($codon{$triplet}) ? $codon{$triplet} : "X";
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
	my $id_order = 0;
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
				"id_order" => $id_order,
				"id_start" => $id_start,
				"seq_length" => $seq_length,
				"seq_start" => $seq_start,
				"row_width" => $row_width,
				"row_bytes" => $row_bytes
			};
			
			# ID開始点を更新
			$id_start = $seq_start + int($seq_length / $row_width) * $row_bytes + $seq_length % $row_width;
			
			# 行番号を加算
			$id_order++;
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
					"id_order" => $id_order,
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
	
	# インデックスハッシュのリファレンスを返す
	return(\%fasta_index);
}

# アイソフォームを決定 common::define_isoforms(bedデータリストリファレンス, 出力アイソフォーム数, 出力遺伝子型)
sub define_isoforms {
	# 変数を宣言
	my %type = ("red" => substr($_[2], 0, 1), "yellow" => substr($_[2], 1, 1), "blue" => substr($_[2], 2, 1), "." => 1);
	my @loci = ();
	
	# データを順に整理
	foreach my $bed (@{$_[0]}) {
		# 未指定の遺伝子タイプのデータを除外
		if (!$type{$bed->[8]}) {next;}
		
		# ブロック情報を絶対位置に変換して最後列に追加
		my @blocks = ();
		for (my $i = 0;$i < $bed->[9];$i++) {$blocks[$i] = [$bed->[1] + $bed->[11]->[$i], $bed->[1] + $bed->[10]->[$i] + $bed->[11]->[$i]];}
		
		# 変数を宣言
		my $basal_locus = "";
		
		# 既出領域と比較
		for (my $i = 0;$i < @loci;$i++) {
			# 方向が異なる場合を除外
			if ($bed->[5] != $loci[$i]->[0]->[3]) {next;}
			
			# 領域がオーバーラップしない場合を除外
			if ($bed->[2] <= $loci[$i]->[0]->[1] or $loci[$i]->[0]->[2] <= $bed->[1]) {next;}
			
			# ブロックがオーバーラップしない場合を除外
			if (!grep {my $block = $_;List::Util::first {$_->[0] < $block->[1] and $block->[0] < $_->[1]} @{$loci[$i]->[0]->[4]}} @blocks) {next;}
			
			# 既存アイソフォームから完全に一致するものを探索
			my $synonym = List::Util::first {join(",", ($bed->[1], @{$bed->[10]}, @{$bed->[11]})) eq join(",", ($_->[1], @{$_->[10]}, @{$_->[11]}))} @{$loci[$i]->[1]};
			
			# 既存アイソフォームと完全に一致する場合を除外
			if ($synonym) {
				if ($bed->[4] > $synonym->[4] or $bed->[4] == $synonym->[4] and $bed->[3] lt $synonym->[3]) {($synonym->[3], $synonym->[4]) = ($bed->[3], $bed->[4]);}
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
	
	# 並べ替えた (開始点 > 終了点 > 向き) 領域リストのリファレンスを返す
	return([sort {$a->[0]->[1] <=> $b->[0]->[1] or $a->[0]->[2] <=> $b->[0]->[2] or $a->[0]->[3] cmp $b->[0]->[3]} @loci]);
}

# データを出力 common::output(bedデータリストリファレンス, 遺伝子番号, プレフィックス, 出力形式)
sub output {
	# 変数を宣言
	my %type = ("red" => "pseudogene", "yellow" => "truncated", "blue" => "protein_coding", "." => "uncharacterized");
	
	# データを出力
	foreach my $locus (@{$_[0]}) {
		# 領域番号を加算
		$_[1]++;
		
		# 先頭の遺伝子行を抽出
		my $gene = shift(@{$locus});
		
		# gtf形式の場合
		if ($_[3] eq "gtf") {
			# 開始点をgtf形式に修正
			$gene->[1]++;
			
			# 遺伝子行を出力
			print "$gene->[0]\tfate\tgene\t$gene->[1]\t$gene->[2]\t.\t", $gene->[3] > 0 ? "+" : "-", "\t.\t", 'gene_id "', "$_[2]$_[1]", '";', "\n";
		}
		
		# アイソフォームを順に処理
		for (my $i = 0;$i < @{$locus};$i++) {
			# bed形式の場合
			if ($_[3] eq "bed") {
				$locus->[$i]->[3] = "$_[2]$_[1]." . ($i + 1) . ":" . substr($locus->[$i]->[3], index($locus->[$i]->[3], ":") + 1);
				$locus->[$i]->[5] = $locus->[$i]->[5] > 0 ? "+" : "-";
				print join("\t", @{$locus->[$i]}[0..9]), "\t", join(",", @{$locus->[$i]->[10]}), "\t", join(",", @{$locus->[$i]->[11]}), "\n";
			}
			
			# gtf形式の場合
			elsif ($_[3] eq "gtf") {
				# 開始点をgtf形式に修正
				$locus->[$i]->[1]++;
				
				# 領域を示す行を出力
				print "$locus->[$i]->[0]\tfate\ttranscript\t$locus->[$i]->[1]\t$locus->[$i]->[2]\t$locus->[$i]->[4]\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $type{$locus->[$i]->[8]}, '";', "\n";
				
				# 逆鎖の場合は個々のブロックを逆順に並べ替え
				if ($locus->[$i]->[5] < 0) {$locus->[$i]->[10] = [reverse(@{$locus->[$i]->[10]})];$locus->[$i]->[11] = [reverse(@{$locus->[$i]->[11]})];}
				
				# 向きを登録
				my $sign = $locus->[$i]->[5];
				$locus->[$i]->[5] = $locus->[$i]->[5] > 0 ? "+" : "-";
				
				# 個々のヒットを示す行を出力
				for (my $j = 0;$j < $locus->[$i]->[9];$j++) {
					# エキソン行を出力
					print "$locus->[$i]->[0]\tfate\texon\t", $locus->[$i]->[1] + $locus->[$i]->[11]->[$j], "\t", $locus->[$i]->[1] + $locus->[$i]->[10]->[$j] + $locus->[$i]->[11]->[$j] - 1, "\t.\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; exon_number "', $j + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $type{$locus->[$i]->[8]}, '"; exon_id "', "$_[2]$_[1].", $i + 1, ".", $j + 1, '";', "\n";
					
					# コーディング領域以外を除外
					if ($type{$locus->[$i]->[8]} ne "protein_coding") {next;}
					
					# コーディング領域行を出力
					print "$locus->[$i]->[0]\tfate\tCDS\t", $locus->[$i]->[1] + $locus->[$i]->[11]->[$j] + 3 * 0 ** ($locus->[$i]->[9] - $j + $sign), "\t", $locus->[$i]->[1] + $locus->[$i]->[10]->[$j] + $locus->[$i]->[11]->[$j] - 1 - 3 * 0 ** ($locus->[$i]->[9] - $j - $sign), "\t.\t$locus->[$i]->[5]\t.\t", 'gene_id "', "$_[2]$_[1]", '"; transcript_id "', "$_[2]$_[1].", $i + 1, '"; exon_number "', $j + 1, '"; transcript_name "', $locus->[$i]->[3], '"; transcript_biotype "', $type{$locus->[$i]->[8]}, '"; protein_id "', "$_[2]$_[1].", $i + 1, 'p";', "\n";
				}
			}
		}
	}
	return(1);
}

# サブルーチンを追加

# パッケージを追加
### 編集範囲 終了 ###
