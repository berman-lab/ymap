<?php
function process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput) {
fwrite($logOutput, "\tPHP : Process uploaded data files into standard forms for pipeline use.\n");
fwrite($logOutput, "\t\t*========================================================*\n");
fwrite($logOutput, "\t\t| Log of 'process_input_files.php'                       |\n");
fwrite($logOutput, "\t\t*--------------------------------------------------------*\n");
fwrite($logOutput, "\t\t| Before archive decompression.\n");
fwrite($logOutput, "\t\t|\text         = ".$ext."\n");
fwrite($logOutput, "\t\t|\tname        = ".$name."\n");
fwrite($logOutput, "\t\t|\tprojectPath = ".$projectPath."\n");
// Deal with compressed archives.
if (strcmp($ext,"zip") == 0) {
	fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
	fwrite($logOutput, "\t\t| This is a ZIP archive of : \n");
	$currentDir = getcwd();                             // get script's path.
	// figure out filename contained in zip archive.
	$null               = shell_exec("unzip -l ".$projectPath.$name." > ".$projectPath."zipTemp.txt");   // generate txt file containing archive contents.
	$zipTempLines       = file($projectPath."zipTemp.txt");
	$zipTempArchiveLine = trim($zipTempLines[3]);
	$columns            = preg_split('/\s+/', $zipTempArchiveLine);
	$oldName            = $columns[3];
	$fileName_parts     = preg_split('/[.]/', $oldName);
	fwrite($logOutput,"\t\t|\t'".$oldName."'.\n");
	fwrite($logOutput,"\t\t| The number of fileName parts = ".count($fileName_parts)."\n");
	// extract archive.
	chdir($projectPath);                   // move to projectDirectory.
	$null = shell_exec("unzip -j ".$name); // unzip archive.
	chdir($currentDir);                    // move back to script's path.
	// delete original archive.
	unlink($projectPath.$name);
	// figure out file type.
	if (count($fileName_parts) == 2) {
		$ext_new        = end($fileName_parts);
	} else {
		// A file extension is not found.
		// determine if it is a FASTQ file by looking at first four lines of text.
		$currentDir = getcwd();
		//fwrite($logOutput, "\t\t| Current directory = '$currentDir'.\n");
		$file_name = $currentDir."/".$projectPath.$oldName;
		//fwrite($logOutput, "\t\t| Open file '".$file_name."'.\n");
		$file_handle = fopen($file_name,'r');
		//fwrite($logOutput, "\t\t| file_handle = '".$file_handle."'\n");
		$line_1      = fgets($file_handle);
		$line_2      = fgets($file_handle);
		$line_3      = fgets($file_handle);
		$line_4      = fgets($file_handle);
		fclose($file_handle);
		if (($line_1[0] == '@') && ($line_3[0] == '+')) {
			// is a FASTQ file.
			$ext_new = "fastq";
		} else {
			// Determine if the file is TXT or array(XLS) by looking at first line.
			$line_1_words = explode(' ',trim($line_1));
			if ((strcmp($line_1_words[0],"Ca19-mtDNA") == 0) || (strcmp($line_1_words[0],"Ca21chr1_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr2_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr3_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr4_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr5_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr6_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr7_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chrR_C_albicans_SC5314") == 0)) {
				// is a TXT file, as defined in the paper.
				$ext_new = "txt";
			} elseif (strcmp($line_1_words[0],"Created") == 0) {
				// is a XLS or a Tab-delimited-txt (TXT) file, as output from BlueFuse.
				$ext_new = "xls";
			} else {
				$ext_new = "none1";
			}
		}
	}
	$name_new           = $oldName;
	// rename decompressed file.
	$rename      = "datafile_".$key.".".$ext_new;
	chdir($projectPath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new    = $rename;
} else if (strcmp($ext,"gz") == 0) {
	fwrite($condensedLogOutput, "Decompressing GZ file : ".$name_new."\n");
	fwrite($logOutput, "\t\t| This is a GZ archive of : \n");
	// figure out filename contained in zip archive.
	$null               = shell_exec("gzip -l ".$projectPath.$name." > ".$projectPath."gzTemp.txt");   // generate txt file containing archive contents.
	$gzTempLines        = file($projectPath."gzTemp.txt");
	$gzTempArchiveLine  = trim($gzTempLines[1]);
	$columns            = preg_split('/\s+/', $gzTempArchiveLine);
	$oldName            = str_replace($projectPath,'',$columns[3]);  // clean off directory path from calling script to datafile.
	$fileName_parts     = preg_split('/[.]/', $oldName);
	fwrite($logOutput,"\t\t|\t'".$oldName."'.\n");
	// extract archive.
	$currentDir = getcwd();                // get script's path.
	chdir($projectPath);                   // move to projectDirectory.
	$null = shell_exec("gzip -d ".$name);  // decompress archive.
	chdir($currentDir);                    // move back to script's path.
	// Figure out file type.
	$currentDir      = getcwd();
	fwrite($logOutput, "\t\t| Current directory  = '".$currentDir."'.\n");
	$file_name       = $currentDir."/".$projectPath.$oldName;
	fwrite($logOutput, "\t\t| Open file '".$file_name."'.\n");
	if (count($fileName_parts) >= 2) {
		// A file extension is found.
		// For example:
		//    [0] = 'test';
		//    [1] = 'fastq';
		$ext_new         = end($fileName_parts);
		fwrite($logOutput, "\t\t| Original file name = '".$oldName."'.\n");
		fwrite($logOutput, "\t\t|    File extension  = '".$ext_new."'\n");
	} else {
		// A file extension is not found.
		// For example:
		//    [0] = 'test';
		// determine if it is a FASTQ file by looking at first four lines of text.
		$file_handle = fopen($file_name,'r');
		$line_1      = fgets($file_handle);
		$line_2      = fgets($file_handle);
		$line_3      = fgets($file_handle);
		$line_4      = fgets($file_handle);
		fclose($file_handle);
		if (($line_1[0] == '@') && ($line_3[0] == '+')) {
			// is a FASTQ file.
			$ext_new = "fastq";
		} else {
			// Determine if the file is TXT or array(XLS) by looking at first line.
			$line_1_words = explode(' ',trim($line_1));
			if ((strcmp($line_1_words[0],"Ca19-mtDNA") == 0) || (strcmp($line_1_words[0],"Ca21chr1_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr2_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr3_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr4_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr5_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chr6_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr7_C_albicans_SC5314") == 0)
				|| (strcmp($line_1_words[0],"Ca21chrR_C_albicans_SC5314") == 0)) {
				// is a TXT file, as defined in the paper. This will only catch datasets for Candida albicans Assembly 21.
				$ext_new = "txt";
			} elseif (strcmp($line_1_words[0],"Created") == 0) {
				// is a XLS or a Tab-delimited-txt (TXT) file, as output from BlueFuse.
				$ext_new = "xls";
			} else {
				$ext_new = "none1";
				// unable to determine file type.
			}
		}
		fwrite($logOutput, "\t\t| Original file name         = '".$oldName."'.\n");
		fwrite($logOutput, "\t\t|    Inferred file extension = '".$ext_new."'\n");
		$oldName = $oldName.".".$ext_new;
		fwrite($logOutput, "\t\t|    New file name           = '".$oldName."'\n");
	}
	// rename decompressed file.
	$rename      = "datafile_".$key.".".$ext_new;
	chdir($projectPath);
	fwrite($logOutput, "\t\t| oldName = '".$oldName."'\n");
	fwrite($logOutput, "\t\t| rename  = '".$rename."'\n");
	rename($oldName,$rename);
	chdir($currentDir);
	$name_new    = $rename;
} else {
	// extension indicates the file is not a compressed archive.
	$ext_new  = $ext;
	$name_new = $name;
}
if (strcmp($ext_new,"") == 0) {
	// A file extension is not found.
	// determine if it is a FASTQ file by looking at first four lines of text.
	$currentDir = getcwd();
	fwrite($logOutput, "\t\t| Current directory = '$currentDir'.\n");
	$file_name = $currentDir."/".$projectPath.$name;
	fwrite($logOutput, "\t\t| Open file '".$file_name."'.\n");
	$file_handle = fopen($file_name,'r');
	fwrite($logOutput, "\t\t| file_handle = '".$file_handle."'\n");
	$line_1      = fgets($file_handle);
	$line_2      = fgets($file_handle);
	$line_3      = fgets($file_handle);
	$line_4      = fgets($file_handle);
	fclose($file_handle);
	if (($line_1[0] == '@') && ($line_3[0] == '+')) {
		// is a FASTQ file.
		$ext_new = "fastq";
	} else {
		// Determine if the file is TXT or array(XLS) by looking at first line.
		$line_1_words = explode(' ',trim($line_1));
		if ((strcmp($line_1_words[0],"Ca19-mtDNA") == 0) || (strcmp($line_1_words[0],"Ca21chr1_C_albicans_SC5314") == 0)
		|| (strcmp($line_1_words[0],"Ca21chr2_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr3_C_albicans_SC5314") == 0)
		|| (strcmp($line_1_words[0],"Ca21chr4_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr5_C_albicans_SC5314") == 0)
		|| (strcmp($line_1_words[0],"Ca21chr6_C_albicans_SC5314") == 0) || (strcmp($line_1_words[0],"Ca21chr7_C_albicans_SC5314") == 0)
		|| (strcmp($line_1_words[0],"Ca21chrR_C_albicans_SC5314") == 0)) {
			// is a TXT file, as defined in the paper.
			$ext_new = "txt";
		} elseif (strcmp($line_1_words[0],"Created") == 0) {
			// is a XLS or a Tab-delimited-txt (TXT) file, as output from BlueFuse.
			$ext_new = "xls";
		} else {
			$ext_new = "none2";
		}
	}
	$name_new    = $name;
	// rename uploaded file with no extension.
	$rename      = "datafile_".$key.".".$ext_new;
	chdir($projectPath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new    = $rename;
}
fwrite($logOutput, "\t\t| After archive decompression.\n");
fwrite($logOutput, "\t\t|\text_new     = ".$ext_new."\n");
fwrite($logOutput, "\t\t|\tname_new    = ".$name_new."\n");
fwrite($logOutput, "\t\t|\tprojectPath = ".$projectPath."\n");

// Deal with raw/decompressed files.
if ((strcmp($ext_new,"fastq") == 0) || (strcmp($ext_new,"fq") == 0)) {
	fwrite($logOutput, "\t\t| This is an uncompressed FASTQ file, no further pre-processing is needed.\n");
	fwrite($output, $name_new."\n");
	$paired = 0;
} elseif ((strcmp($ext_new,"fasta") == 0) || (strcmp($ext_new,"fa") == 0)) {
	fwrite($logOutput, "\t\t| This is a FASTA file.\n");
	$errorFile = fopen("users/".$user."/projects/".$project."/error.txt", 'w');
	fwrite($errorFile, "Error : FASTA file uploaded as input. Upload FASTQ, or ZIP or GZ archives.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
} elseif ((strcmp($ext_new,"sam") == 0) || (strcmp($ext_new,"bam") == 0)) {
	fwrite($logOutput, "\t\t| This is a SAM/BAM file.\n");

	// The .sh scripts assume we're running from Ymap root.
    $absProjectPath = realpath($projectPath) . "/";
	$currentDir = getcwd();
	chdir($projectPath . "/../../../../"); // ymap_root/users/user_name/projects/this_project

	// Convert BAM to SAM file, if needed.
	if (strcmp($ext_new,"bam") == 0) {
		fwrite($condensedLogOutput, "Decompressing BAM file to SAM.\n");
		$null = shell_exec("sh scripts_seqModules/bam2sam.sh ".$user." ".$project." ".$name_new);
		unlink($absProjectPath.$name_new);
		unlink($absProjectPath.$name_new.".bai");
		$name_new = "data.sam";
	}
	// Convert SAM file to FASTQ files.
	fwrite($condensedLogOutput, "Decompressing SAM file to FASTQ.\n");
	$null       = shell_exec("sh scripts_seqModules/sam2fastq.sh ".$user." ".$project." ".$name_new);
	// sam2fastq.sh user project main_dir inputFile;
	fwrite($output, "data_r1.fastq\n");
	fwrite($output, "data_r2.fastq\n");
	// delete original archive.
	unlink($absProjectPath.$name_new);
	fwrite($logOutput, "\t\t| File converted to paired-FASTQ files, original deleted.\n");
	$paired = 1;

	chdir($currentDir);
} elseif (strcmp($ext_new,"txt") == 0) {
	fwrite($logOutput, "\t\t| This is a txt file.\n");
	fwrite($logOutput, "\t\t|\tCurrentDir = ".getcwd()."\n");
	fwrite($logOutput, "\t\t|\tshell_exec string = 'sh ../Gareth2pileups.sh ".$user." ".$project." ".$name_new."'\n");
	$currentDir = getcwd();
	$null       = shell_exec("sh ../Gareth2pileups.sh ".$user." ".$project." ".$name_new);
	// sam2fastq.sh user project main_dir inputFile;
	fwrite($output, "null1\n");
	fwrite($output, "null2\n");
	// delete original archive.
	unlink($projectPath.$name_new);
	fwrite($logOutput, "\t\t| File converted to intermediate pileup files, original deleted.\n");
	$paired = 1;
} elseif (strcmp($ext_new,"none1") == 0) {
	fwrite($logOutput, "\t\t| This archive contained a file with no extension and the file type could not be determined.\n");
	$errorFile = fopen("users/".$user."/projects/".$project."/error.txt", 'w');
	fwrite($errorFile, "Error : Archive contained a file with no extension and the file type could not be determined.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
} elseif (strcmp($ext_new,"none2") == 0) {
	fwrite($logOutput, "\t\t| This file had no extension and the file type could not be determined.\n");
	$errorFile = fopen("users/".$user."/projects/".$project."/error.txt", 'w');
	fwrite($errorFile, "Error : File had no extension and the file type could not be determined.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
} else {
	fwrite($logOutput, "\t\t| This is an unknown file type.\n");
	$errorFile = fopen("users/".$user."/projects/".$project."/error.txt", 'w');
	fwrite($errorFile, "Error : Unknown file type as input.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
}

fwrite($logOutput, "\t\t*--------------------------------------------------------*\n");
fwrite($logOutput, "\t\t| 'process_input_files.php' has completed.               |\n");
fwrite($logOutput, "\t\t*========================================================*\n");
return $paired;
}
?>
