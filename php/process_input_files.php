<?php
function process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput) {
fwrite($logOutput, "Before archive decompression.\n");
fwrite($logOutput, "\t| ext         = ".$ext."\n");
fwrite($logOutput, "\t| name        = ".$name."\n");
fwrite($logOutput, "\t| projectPath = ".$projectPath."\n");
	// Deal with compressed archives.
	if (strcmp($ext,"zip") == 0) {
		fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
		fwrite($logOutput, "\t| This is a ZIP archive of : ");
		$currentDir = getcwd();                             // get script's path.
		// figure out filename contained in zip archive.
		$null               = shell_exec("unzip -l ".$projectPath.$name." > ".$projectPath."zipTemp.txt");   // generate txt file containing archive contents.
		$zipTempLines       = file($projectPath."zipTemp.txt");
		$zipTempArchiveLine = trim($zipTempLines[3]);
		$columns            = preg_split('/\s+/', $zipTempArchiveLine);
		$oldName            = $columns[3];
		fwrite($logOutput, " '$oldName'.\n");
		$fileName_parts     = preg_split('/[.]/', $oldName);
		fwrite($logOutput,"\t| number of fileName parts = ".count($fileName_parts)."\n");
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
			$current_Dir = getcwd();
			fwrite($logOutput, "Current directory = '$current_Dir'.\n");
			$file_name = $current_Dir."/".$projectPath.$oldName;
			fwrite($logOutput, "Open file '".$file_name."'.\n");
			$file_handle = fopen($file_name,'r');
			fwrite($logOutput, "file_handle = '".$file_handle."'\n");
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
		fwrite($logOutput, "\t| This is a GZ archive of : ");
		// figure out filename contained in zip archive.
		$null               = shell_exec("gzip -l ".$projectPath.$name." > ".$projectPath."gzTemp.txt");   // generate txt file containing archive contents.
		$gzTempLines        = file($projectPath."gzTemp.txt");
		$gzTempArchiveLine  = trim($gzTempLines[1]);
		$columns            = preg_split('/\s+/', $gzTempArchiveLine);
		$oldName            = $columns[3];
		$fileName_parts     = preg_split('/[.]/', $oldName);
		fwrite($logOutput, " '$oldName'.\n");
		fwrite($logOutput,"\t| number of fileName parts = ".count($fileName_parts)."\n");
		// extract archive.
		$currentDir = getcwd();                // get script's path.
		chdir($projectPath);                   // move to projectDirectory.
		$null = shell_exec("gzip -d ".$name);  // decompress archive.
		chdir($currentDir);                    // move back to script's path.
		// Figure out file type.
		if (count($fileName_parts) == 4) {
			// A file extension is found.
			// For example:
			//    [0] = '.';
			//    [1] = '.';
			//    [2] = '/users/darren1/projects/test_gz4/test';
			//    [3] = 'fastq';
			$ext_new        = end($fileName_parts);
		} else {
			// A file extension is not found.
			// determine if it is a FASTQ file by looking at first four lines of text.
			$current_Dir = getcwd();
            fwrite($logOutput, "Current directory = '$current_Dir'.\n");
			$file_name = $current_Dir."/".$oldName;
			fwrite($logOutput, "Open file '".$file_name."'.\n");
			$file_handle = fopen($file_name,'r');
			fwrite($logOutput, "file_handle = '".$file_handle."'\n");
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
	} else {
		$ext_new  = $ext;
		$name_new = $name;
	}
if (strcmp($ext_new,"") == 0) {
	// A file extension is not found.
	// determine if it is a FASTQ file by looking at first four lines of text.
	$current_Dir = getcwd();
	fwrite($logOutput, "Current directory = '$current_Dir'.\n");
	$file_name = $current_Dir."/".$projectPath.$name;
	fwrite($logOutput, "Open file '".$file_name."'.\n");
	$file_handle = fopen($file_name,'r');
	fwrite($logOutput, "file_handle = '".$file_handle."'\n");
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
}
fwrite($logOutput, "After archive decompression.\n");
fwrite($logOutput, "\t| ext_new     = ".$ext_new."\n");
fwrite($logOutput, "\t| name_new    = ".$name_new."\n");
fwrite($logOutput, "\t| projectPath = ".$projectPath."\n");

	// Deal with raw/decompressed files.
	if ((strcmp($ext_new,"fastq") == 0) || (strcmp($ext_new,"fq") == 0)) {
		fwrite($logOutput, "\t| This is an uncompressed FASTQ file, no further pre-processing is needed.\n");
		fwrite($output, $name_new."\n");
	} elseif ((strcmp($ext_new,"fasta") == 0) || (strcmp($ext_new,"fa") == 0)) {
		fwrite($logOutput, "\t| This is a FASTA file.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : FASTA file uploaded as input. Upload FASTQ, or ZIP or GZ archives.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	} elseif ((strcmp($ext_new,"sam") == 0) || (strcmp($ext_new,"bam") == 0)) {
		fwrite($logOutput, "\t| This is a SAM/BAM file.\n");
		// Convert BAM to SAM file, if needed.
		if (strcmp($ext_new,"bam") == 0) {
			fwrite($condensedLogOutput, "Decompressing BAM file to SAM.\n");
			$null = shell_exec("sh ../sh/bam2sam.sh ".$user." ".$project." ".$name_new);
			unlink($projectPath.$name_new);
			unlink($projectPath.$name_new.".bai");
			$name_new = "data.sam";
		}
		// Convert SAM file to FASTQ files.
		fwrite($condensedLogOutput, "Decompressing SAM file to FASTQ.\n");
		$currentDir = getcwd();
		$null       = shell_exec("sh ../sh/sam2fastq.sh ".$user." ".$project." ".$name_new);
		// sam2fastq.sh user project main_dir inputFile;
		fwrite($output, "data_r1.fastq\n");
		fwrite($output, "data_r2.fastq\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to FASTQ files, original deleted.\n");
		$paired = 1;
	} elseif (strcmp($ext_new,"txt") == 0) {
		fwrite($logOutput, "\t| This is a txt file.\n");
		$currentDir = getcwd();
		$null       = shell_exec("sh ../sh/Gareth2pileups.sh ".$user." ".$project." ".$name_new);
		// sam2fastq.sh user project main_dir inputFile;
		fwrite($output, "null1\n");
		fwrite($output, "null2\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to FASTQ files, original deleted.\n");
		$paired = 1;
	} elseif (strcmp($ext_new,"none1") == 0) {
		fwrite($logOutput, "\t| This archive contained a file with no extension and the file type could not be determined.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : Archive contained a file with no extension and the file type could not be determined.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	} elseif (strcmp($ext_new,"none2") == 0) {
		fwrite($logOutput, "\t| This file had no extension and the file type could not be determined.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : File had no extension and the file type could not be determined.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	} else {
		fwrite($logOutput, "\t| This is an unknown file type.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : Unknown file type as input.\nUpload FASTQ, or ZIP or GZ archives containing a FASTQ file.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	}

	return $paired;
}
?>
