<?php
function process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput) {
fwrite($logOutput, "Before archive decompression.\n");
fwrite($logOutput, "\t| ext         = ".$ext."\n");
fwrite($logOutput, "\t| name        = ".$name."\n");
fwrite($logOutput, "\t| projectPath = ".$projectPath."\n");
	// Deal with compressed archives.
	if (strcmp($ext,"zip") == 0) {
		fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
		fwrite($logOutput, "\t\tThis is a ZIP archive of : ");
		$currentDir = getcwd();                             // get script's path.
		// figure out filename contained in zip archive.
		$null               = shell_exec("unzip -l ".$projectPath.$name." > ".$projectPath."zipTemp.txt");   // generate txt file containing archive contents.
		$zipTempLines       = file($projectPath."zipTemp.txt");
		$zipTempArchiveLine = trim($zipTempLines[3]);
		$columns            = preg_split('/\s+/', $zipTempArchiveLine);
		$oldName            = $columns[3];
		fwrite($logOutput, " '$oldName'.\n");
		$fileName_parts     = preg_split('/[.]/', $oldName);
		$ext_new            = end($fileName_parts);
		$name_new           = $oldName;
		// extract archive.
		chdir($projectPath);                   // move to projectDirectory.
		$null = shell_exec("unzip -j ".$name); // unzip archive.
		chdir($currentDir);                    // move back to script's path.
		// delete original archive.
		unlink($projectPath.$name);
		// rename decompressed file.
		$rename      = "datafile_".$key.".".$ext_new;
		chdir($projectPath);
		rename($name_new,$rename);
		chdir($currentDir);
		$name_new    = $rename;
	} else if (strcmp($ext,"gz") == 0) {
		fwrite($condensedLogOutput, "Decompressing GZ file : ".$name_new."\n");
		fwrite($logOutput, "\t\tThis is a GZ archive of : ");
		// figure out filename contained in zip archive.
		$null               = shell_exec("gzip -l ".$projectPath.$name." > ".$projectPath."gzTemp.txt");   // generate txt file containing archive contents.
		$gzTempLines        = file($projectPath."gzTemp.txt");
		$gzTempArchiveLine  = trim($gzTempLines[1]);
		$columns            = preg_split('/\s+/', $gzTempArchiveLine);
		$oldName            = $columns[3];
		fwrite($logOutput, " '$oldName'.\n");
		$fileName_parts     = preg_split('/[.]/', $oldName);
		$ext_new            = end($fileName_parts);
		$name_new           = $oldName;
		// extract archive.
		$currentDir = getcwd();                // get script's path.
		chdir($projectPath);                   // move to projectDirectory.
		$null = shell_exec("gzip -d ".$name);  // decompress archive.
		chdir($currentDir);                    // move back to script's path.
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
fwrite($logOutput, "After archive decompression.\n");
fwrite($logOutput, "\t| ext_new     = ".$ext_new."\n");
fwrite($logOutput, "\t| name_new    = ".$name_new."\n");
fwrite($logOutput, "\t| projectPath = ".$projectPath."\n");

	// Deal with raw/decompressed files.
	if ((strcmp($ext_new,"fastq") == 0) || (strcmp($ext_new,"fq") == 0)) {
		fwrite($logOutput, "\t\tThis is an uncompressed FASTQ file, no further pre-processing is needed.\n");
		fwrite($output, $name_new."\n");
	} else if ((strcmp($ext_new,"fasta") == 0) || (strcmp($ext_new,"fa") == 0)) {
		fwrite($logOutput, "\tThis is a FASTA file.\n");
		$errorFile = fopen("../users/".$user."/projects/".$project."/error.txt", 'w');
		fwrite($errorFile, "Error : FASTA file uploaded as input. Upload FASTQ, or ZIP or GZ archives.");
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
		$currentDir = getcwd();
		$null       = shell_exec("sh ../sh/sam2fastq.sh ".$user." ".$project." ".$name_new);
		// sam2fastq.sh user project main_dir inputFile;
		fwrite($output, "data_r1.fastq\n");
		fwrite($output, "data_r2.fastq\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to FASTQ files, original deleted.\n");
		$paired = 1;
	} else if (strcmp($ext_new,"txt") == 0) {
		fwrite($logOutput, "\tThis is a txt file.\n");
		$currentDir = getcwd();
		$null       = shell_exec("sh ../sh/Gareth2pileups.sh ".$user." ".$project." ".$name_new);
		// sam2fastq.sh user project main_dir inputFile;
		fwrite($output, "null1\n");
		fwrite($output, "null2\n");
		// delete original archive.
		unlink($projectPath.$name_new);
		fwrite($logOutput, "\t\tFile converted to FASTQ files, original deleted.\n");
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
