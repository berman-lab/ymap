<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<HTML>
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
		.upload {
			width:          675px;      // 675px;
			border:         0;
			height:         40px;   // 40px;
			vertical-align: middle;
			align:          left;
			margin:         0px;
			overflow:       hidden;
		}
		html, body {
			margin:         0px;
			border:         0;
			overflow:       hidden;
		}
	</style>
<meta http-equiv="content-type" content="text/html; charset=iso-8859-1">
<title>Install project into pipeline.</title>
</HEAD>
<?php
    require_once 'constants.php';

// Deal with passed variables.
	$fileName         = $argv[1];
	$user             = $argv[2];
	$project          = $argv[3];

//$fileName = "Cg_r1.zip,Cg_r2.zip";
//$user     = "darren";
//$project  = "test_Cg";

// Initialize log file.
	$logOutputName = $directory."users/".$user."/projects/".$project."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "#..............................................................................\n");
	fwrite($logOutput, "Running 'php/project.paired_WGseq.install_2.php'.\n");
	fwrite($logOutput, "Variables passed via command-line from 'php/project.paired_WGseq.install_1.php' :\n");
	fwrite($logOutput, "\tfileName         = '".$fileName."'\n");
	fwrite($logOutput, "\tuser             = '".$user."'\n");
	fwrite($logOutput, "\tproject          = '".$project."'\n");
	fwrite($logOutput, "#============================================================================== 3\n");

// Manage condensed log file.
	$condensedLogOutputName = $directory."users/".$user."/projects/".$project."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'a');
//	fclose($condensedLogOutput);

// Generate 'datafiles.txt' file containing: name of all data files.
// Identify format of uploaded file and decompress as needed (*.ZIP; *.GZ).
	$outputName = $directory."users/".$user."/projects/".$project."/datafiles.txt";
	$output     = fopen($outputName, 'w');
	$fileNames = explode(",", $fileName);
	fwrite($logOutput, "\tGenerate 'datafiles.txt' and decompress uploaded archives.\n");
	foreach ($fileNames as $key=>$name) {
		$projectPath = $directory."users/".$user."/projects/".$project."/";
		$name        = str_replace("\\", ",", $name);
		$ext         = strtolower(pathinfo($name, PATHINFO_EXTENSION));
		$filename    = strtolower(pathinfo($name, PATHINFO_FILENAME));
		fwrite($logOutput, "\tFile ".$key."\n");
		fwrite($logOutput, "\t\tDatafile   : '$name'.\n");
		fwrite($logOutput, "\t\tFilename   : '$filename.'.\n");
		fwrite($logOutput, "\t\tExtension  : '$ext'.\n");
		fwrite($logOutput, "\t\tPath       : '$projectPath'.\n");
		// Standardize name of uploaded file.
		$oldName     = $name;
		$newName     = "datafile_".$key.".".$ext;
		chdir($projectPath);
		rename($oldName, $newName);
		chdir($currentPath);
		unlink($projectPath.$name);

		if (strcmp($ext,"zip") == 0) {
			fwrite($condensedLogOutput, "Decompressing ZIP file : ".$name."\n");
			fwrite($logOutput, "\t\tThis is a ZIP file.\n");
			$currentDir = getcwd(); 							// get script's path.
			chdir($projectPath);								// move to projectDirectory.
			$null = shell_exec("unzip -j ".$projectPath.$name);	// unzip archive.
			chdir($currentDir);									// move back to script's path.
			// figure out filename contained in zip archive.
			$null               = shell_exec("unzip -l ".$projectPath.$name." > ".$projectPath."zipTemp.txt");   // generate txt file containing archive contents.
			$zipTempLines       = file($projectPath."zipTemp.txt");
			$zipTempArchiveLine = $zipTempLines[3];
			$columns            = preg_split('/\s+/', $zipTempArchiveLine);
			$oldName            = $columns[3];
			$newName            = "datafile_".$key.".fastq";
			chdir($projectPath);
			rename($oldName, $newName);
			chdir($currentDir);
			fwrite($logOutput, "\t\tArchive of : '$newName'.\n");
			// delete original archive.
			unlink($projectPath.$name);
			fwrite($logOutput, "\t\tFile unzipped, original deleted.\n");
			// add file name to 'datafiles.txt'.
			fwrite($output, $newName."\n");
		} else if (strcmp($ext,"gz") == 0) {
			fwrite($condensedLogOutput, "Decompressing GZ file : ".$name."\n");
			fwrite($logOutput, "\t\tThis is a GZ file.\n");
			$currentDir = getcwd(); 							// get script's path.
			chdir($projectPath);								// move to projectDirectory.
			$null = shell_exec("gzip -d ".$projectPath.$name);	// decompress archive.
			chdir($currentDir);									// move back to script's path.
			$oldName = str_replace(".gz","",$name);
			$newName = "datafile_".$key.".fastq";
			chdir($projectPath);
			rename($oldName,$newName);
			chdir($currentDir);
			fwrite($logOutput, "\t\tFile uncompressed, original deleted.\n");
			// add file name to 'datafiles.txt'.
			fwrite($output, $newName."\n");
		} else if (strcmp($ext,"fastq") == 0) {
			fwrite($logOutput, "\t\tThis is an uncompressed FASTQ file.\n");
			fwrite($output, $name."\n");
		} else {
			fwrite($logOutput, "\tThis is an unknown file type.\n");

			$errorFile = fopen($directory."users/".$user."/projects/".$project."/error.txt", 'w');
			fwrite($errorFile, "Error : Uknown file type as input. Upload FASTQ, or ZIP or GZ archives.");
			fclose($errorFile);
			chmod($errorFileName,0755);
			exit;
		}
	}
	fclose($output);
	chmod($outputName,0755);
	// Trim the last "\n" character from 'datafiles.txt'.
	$fh = fopen($outputName, 'r+');
	$stat = fstat($fh);
	ftruncate($fh, $stat['size']-1);
	fclose($fh);
	fwrite($logOutput, "Completed 'datafiles.txt' file.\n");

	// Final install functions are in shell script.
	fwrite($logOutput, "Passing control to : 'sh/project.paired_WGseq.install_3.sh'\n");
	fwrite($logOutput, "Current directory = '".getcwd()."'\n" );
	$system_call_string = "sh ".$directory."sh/project.paired_WGseq.install_3.sh ".$user." ".$project." ".$directory." > /dev/null &";
	system($system_call_string);
	fclose($condensedLogOutput);
	fclose($logOutput);
?>
