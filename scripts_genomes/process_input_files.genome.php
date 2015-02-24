<?php
function process_input_files_genome($ext,$name,$genomePath,$key,$user,$genome,$output, $condensedLogOutput,$logOutput, $output_fileName) {
fwrite($logOutput, "Before archive decompression.\n");
fwrite($logOutput, "\t| ext         = ".$ext."\n");
fwrite($logOutput, "\t| name        = ".$name."\n");
fwrite($logOutput, "\t| genomePath  = ".$genomePath."\n");
// Deal with compressed archives.
if (strcmp($ext,"zip") == 0) {
	fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
	fwrite($logOutput, "\t| This is a ZIP archive of : ");
	$currentDir = getcwd();                             // get script's path.
	// figure out filename contained in zip archive.
	$null               = shell_exec("unzip -l ".$genomePath.$name." > ".$genomePath."zipTemp.txt");   // generate txt file containing archive contents.
	$zipTempLines       = file($genomePath."zipTemp.txt");
	$zipTempArchiveLine = trim($zipTempLines[3]);
	$columns            = preg_split('/\s+/', $zipTempArchiveLine);
	$oldName            = $columns[3];
	fwrite($logOutput, " '$oldName'.\n");
	$fileName_parts     = preg_split('/[.]/', $oldName);
	fwrite($logOutput,"\t| number of fileName parts = ".count($fileName_parts)."\n");
	// extract archive.
	chdir($genomePath);                   // move to genomeDirectory.
	$null = shell_exec("unzip -j ".$name); // unzip archive.
	chdir($currentDir);                    // move back to script's path.
	// delete original archive.
	unlink($genomePath.$name);
	// figure out file type.
	if (count($fileName_parts) == 2) {
		$ext_new        = end($fileName_parts);
	} else {
		// A file extension is not found.
		// determine if it is a FASTA file by looking at first four lines of text.
		$currentDir = getcwd();
		//fwrite($logOutput, "Current directory = '$currentDir'.\n");
		$file_name = $currentDir."/".$genomePath.$oldName;
		//fwrite($logOutput, "Open file '".$file_name."'.\n");
		$file_handle = fopen($file_name,'r');
		//fwrite($logOutput, "file_handle = '".$file_handle."'\n");
		$line_1      = fgets($file_handle);
		$line_2      = fgets($file_handle);
		$line_3      = fgets($file_handle);
		$line_4      = fgets($file_handle);
		fclose($file_handle);
		if ($line_1[0] == '>') {
			// is a FASTA file.  This is a very poor criteria, but will fail many other formats.
			$ext_new = "fasta";
		} else {
			$ext_new = "none1";
		}
	}
	$name_new           = $oldName;
	// rename decompressed file.
	$rename      = $output_fileName;  //"datafile_".$key.".".$ext_new;
	chdir($genomePath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new    = $rename;
} else if (strcmp($ext,"gz") == 0) {
	fwrite($condensedLogOutput, "Decompressing GZ file : ".$name_new."\n");
	fwrite($logOutput, "\t| This is a GZ archive of : ");
	// figure out filename contained in gz archive.
	$null               = shell_exec("gzip -l ".$genomePath.$name." > ".$genomePath."gzTemp.txt");   // generate txt file containing archive contents.
	$gzTempLines        = file($genomePath."gzTemp.txt");
	$gzTempArchiveLine  = trim($gzTempLines[1]);
	$columns            = preg_split('/\s+/', $gzTempArchiveLine);
	$oldName            = $columns[3];
	$fileName_parts     = preg_split('/[.]/', $oldName);
	fwrite($logOutput, " '$oldName'.\n");
	fwrite($logOutput,"\t| number of fileName parts = ".count($fileName_parts)."\n");
	// extract archive.
	$currentDir = getcwd();                // get script's path.
	chdir($genomePath);                    // move to genomeDirectory.
	$null = shell_exec("gzip -d ".$name);  // decompress archive.
	chdir($currentDir);                    // move back to script's path.
	// Figure out file type.
	if (count($fileName_parts) == 4) {
		// A file extension is found.
		// For example:
		//    [0] = '.';
		//    [1] = '.';
		//    [2] = '/users/darren1/genomes/test_gz4/test';
		//    [3] = 'fastq';
		$fileName_parts2 = preg_split('/[\/]/', $fileName_parts[count($fileName_parts)-2] );
		$oldName         = end($fileName_parts2).".".end($fileName_parts);
		$ext_new         = end($fileName_parts);
	} else {
		// A file extension is not found.
		// For example:
		//    [0] = '.';
		//    [1] = '.';
		//    [2] = '/users/darren1/genomes/test_gz4/test';
		$fileName_parts2 = preg_split('/[\/]/', end($fileName_parts) );
		$oldName         = end($fileName_parts2);
		// determine if it is a FASTQ file by looking at first four lines of text.
		$currentDir      = getcwd();
		// fwrite($logOutput, "Current directory = '$currentDir'.\n");
		$file_name       = $currentDir."/".$genomePath.$oldName;
		// fwrite($logOutput, "Open file '".$file_name."'.\n");
		$file_handle     = fopen($file_name,'r');
		// fwrite($logOutput, "file_handle = '".$file_handle."'\n");
		$line_1          = fgets($file_handle);
		$line_2          = fgets($file_handle);
		$line_3          = fgets($file_handle);
		$line_4          = fgets($file_handle);
		fclose($file_handle);
		if ($line_1[0] == '>') {
			// is a FASTA file.  This is a very poor criteria, but will fail many other formats.
			$ext_new = "fasta";
		} else {
			$ext_new = "none1";
		}
	}
	fwrite($logOutput, "\t| Original name within archive = '".$oldName."'.\n");
	$name_new    = $oldName;
	// rename decompressed file.
	$rename      = $output_fileName;   // "datafile_".$key.".".$ext_new;
	chdir($genomePath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new    = $rename;
} else {
	// extension indicates the file is not a compressed archive.
	$ext_new  = $ext;
	$name_new = $name;

	$rename = $output_fileName;
	$currentDir = getcwd();
	chdir($genomePath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new = $rename;
}
if (strcmp($ext_new,"") == 0) {
	// A file extension is not found.
	// determine if it is a FASTQ file by looking at first four lines of text.
	$currentDir = getcwd();
	fwrite($logOutput, "Current directory = '$currentDir'.\n");
	$file_name = $currentDir."/".$genomePath.$name;
	fwrite($logOutput, "Open file '".$file_name."'.\n");
	$file_handle = fopen($file_name,'r');
	fwrite($logOutput, "file_handle = '".$file_handle."'\n");
	$line_1      = fgets($file_handle);
	$line_2      = fgets($file_handle);
	$line_3      = fgets($file_handle);
	$line_4      = fgets($file_handle);
	fclose($file_handle);
	if ($line_1[0] == '>') {
		// is a FASTA file.  This is a very poor criteria, but will fail many other formats.
		$ext_new = "fasta";
	} else {
		$ext_new = "none2";
	}
	$name_new    = $name;
	// rename uploaded file with no extension.
	$rename      = $output_fileName;   // "datafile_".$key.".".$ext_new;
	chdir($genomePath);
	rename($name_new,$rename);
	chdir($currentDir);
	$name_new    = $rename;
}
fwrite($logOutput, "After archive decompression.\n");
fwrite($logOutput, "\t| ext_new     = ".$ext_new."\n");
fwrite($logOutput, "\t| name_new    = ".$name_new."\n");
fwrite($logOutput, "\t| genomePath  = ".$genomePath."\n");

// Deal with raw/decompressed files.
if ((strcmp($ext_new,"fasta") == 0) || (strcmp($ext_new,"fa") == 0)) {
	fwrite($logOutput, "\t| This is an uncompressed FASTA file, no further pre-processing is needed.\n");
	fwrite($output, $name_new."\n");
} elseif (strcmp($ext_new,"tab") == 0) {
	fwrite($logOutput, "\t| This is an uncompressed TAB file, no further pre-processing is needed.\m");
	fwrite($output, $name_new."\n");
} elseif (strcmp($ext_new,"none1") == 0) {
	fwrite($logOutput, "\t| This archive contained a file with no extension and the file type could not be determined.\n");
	$errorFile = fopen("../users/".$user."/genomes/".$genome."/error.txt", 'w');
	fwrite($errorFile, "Error : Archive contained a file with no extension and the file type could not be determined.\nUpload FASTA, or ZIP or GZ archives containing a FASTA file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
} elseif (strcmp($ext_new,"none2") == 0) {
	fwrite($logOutput, "\t| This file had no extension and the file type could not be determined.\n");
	$errorFile = fopen("../users/".$user."/genomes/".$genome."/error.txt", 'w');
	fwrite($errorFile, "Error : File had no extension and the file type could not be determined.\nUpload FASTA, or ZIP or GZ archives containing a FASTA file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
} else {
	fwrite($logOutput, "\t| This is an unknown file type.\n");
	$errorFile = fopen("../users/".$user."/genomes/".$genome."/error.txt", 'w');
	fwrite($errorFile, "Error : Unknown file type as input.\nUpload FASTA, or ZIP or GZ archives containing a FASTA file.");
	fclose($errorFile);
	chmod($errorFileName,0755);
	exit;
}

return $name_new;
}
?>
