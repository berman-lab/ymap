<?php
function process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput) {
fwrite($logOutput, "\tBefore archive decompression.\n");
fwrite($logOutput, "\t\t| ext         = ".$ext."\n");
fwrite($logOutput, "\t\t| name        = ".$name."\n");
fwrite($logOutput, "\t\t| projectPath = ".$projectPath."\n");
	// Deal with compressed archives.
	$currentDir = getcwd();
	chdir($projectPath);
	if (strcmp($ext,"zip") == 0) {
		fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
		fwrite($logOutput, "\t\t|\tThis is a ZIP archive of : ");
		// figure out filename contained in zip archive.
		$null               = shell_exec("unzip -l ".$name." > zipTemp.txt");   // generate txt file containing archive contents.
		$zipTempLines       = file("zipTemp.txt");
		$zipTempArchiveLine = trim($zipTempLines[3]);
		$columns            = preg_split('/\s+/', $zipTempArchiveLine);
		$oldName            = $columns[3];
		fwrite($logOutput, " '$oldName'.\n");
		fwrite($logOutput, "\t\t|\t\tzipTemp col[0] : '$columns[0]'\n");
		fwrite($logOutput, "\t\t|\t\tzipTemp col[1] : '$columns[1]'\n");
		fwrite($logOutput, "\t\t|\t\tzipTemp col[2] : '$columns[2]'\n");
		fwrite($logOutput, "\t\t|\t\tzipTemp col[3] : '$columns[3]'\n");
		fwrite($logOutput, "\t\t|\t\tzipTemp col[4] : '$columns[4]'\n");
		$fileName_parts     = preg_split('/[.]/', $oldName);
		$ext_new            = end($fileName_parts);
		$name_new           = $oldName;
		// extract archive.
		$null = shell_exec("unzip -j ".$name); // unzip archive.
		// delete original archive.
		unlink($name);
		// rename decompressed file.
		$rename      = "datafile_".$key.".".$ext_new;
		rename($name_new,$rename);
		$name_new    = $rename;
	} else if (strcmp($ext,"gz") == 0) {
		fwrite($condensedLogOutput, "Decompressing GZ file : ".$name_new."\n");
		fwrite($logOutput, "\t\t|\tThis is a GZ archive of : ");
		// figure out filename contained in zip archive.
		$null               = shell_exec("gzip -l ".$name." > gzTemp.txt");   // generate txt file containing archive contents.
		$gzTempLines        = file("gzTemp.txt");
		$gzTempArchiveLine  = trim($gzTempLines[1]);
		$columns            = preg_split('/\s+/', $gzTempArchiveLine);
		$oldName            = $columns[3];
		fwrite($logOutput, " '$oldName'.\n");
		fwrite($logOutput, "\t\t|\t\tgzTemp col[0] : '$columns[0]'\n");
		fwrite($logOutput, "\t\t|\t\tgzTemp col[1] : '$columns[1]'\n");
		fwrite($logOutput, "\t\t|\t\tgzTemp col[2] : '$columns[2]'\n");
		fwrite($logOutput, "\t\t|\t\tgzTemp col[3] : '$columns[3]'\n");
		fwrite($logOutput, "\t\t|\t\tgzTemp col[4] : '$columns[4]'\n");
		$fileName_parts     = preg_split('/[.]/', $oldName);
		$ext_new            = end($fileName_parts);
		// extract archive.
		$null = shell_exec("gzip -d ".$name);  // decompress archive.
		$newName     = "datafile_".$key.".".$ext_new;
		rename($oldName,$newName);
		$name_new    = $newName;
	} else {
		$ext_new  = $ext;
		$name_new = $name;
		$oldName  = $name;
	}
	chdir($currentDir);
fwrite($logOutput, "\tAfter archive decompression.\n");
fwrite($logOutput, "\t\t| ext_new     = ".$ext_new."\n");
fwrite($logOutput, "\t\t| oldName     = ".$oldName."\n");
fwrite($logOutput, "\t\t| name_new    = ".$name_new."\n");
fwrite($logOutput, "\t\t| projectPath = ".$projectPath."\n\n");
$paired = 0;
	// Deal with raw/decompressed files.
	if ((strcmp($ext_new,"fastq") == 0) || (strcmp($ext_new,"fq") == 0)) {
		fwrite($logOutput, "\t\tThis is an uncompressed FASTQ file, no further pre-processing is needed.\n");
		fwrite($output, $name_new."\n");
	} else if ((strcmp($ext_new,"fasta") == 0) || (strcmp($ext_new,"fa") == 0)) {
		fwrite($logOutput, "\tThis is a FASTA file.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : FASTA file uploaded as input, but cannot be processed. Upload FASTQ instead.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	} else if ((strcmp($ext_new,"sam") == 0) || (strcmp($ext_new,"bam") == 0)) {
		fwrite($logOutput, "\tThis is a SAM/BAM file.\n");
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
		$null       = shell_exec("sh ../sh/sam2fastq.sh ".$user." ".$project." ".$name_new);
		fwrite($output, "data_r1.fastq\n");
		fwrite($output, "data_r2.fastq\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to FASTQ files, original deleted.\n");
		$paired = 1;
	} else if (strcmp($ext_new,"txt") == 0) {
		fwrite($logOutput, "\tThis is a txt file.\n");
		// Convert TXT file to pileup files.
		fwrite($condensedLogOutput, "Decompressing TXT file to pileup files.\n");
		$null       = shell_exec("sh ../sh/Gareth2pileups.sh ".$user." ".$project." ".$name_new);
		fwrite($output, "null1\n");
		fwrite($output, "null2\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to pileup files, original deleted.\n");
		$paired = 1;
	} else {
		fwrite($logOutput, "\tThis is an unknown file type.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : Uknown file type as input. Upload FASTQ, or ZIP or GZ archives.");
		fclose($errorFile);
		chmod($errorFileName,0755);
		exit;
	}

	return $paired;
}
?>
