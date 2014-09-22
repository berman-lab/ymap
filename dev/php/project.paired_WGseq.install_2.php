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
	include_once 'process_input_files.php';

// Deal with passed variables.
	$fileName         = $argv[1];
	$user             = $argv[2];
	$project          = $argv[3];

//$fileName = "Cg_r1.zip,Cg_r2.zip";
//$user     = "darren";
//$project  = "test_Cg";

// Initialize log file.
	$logOutputName = "../users/".$user."/projects/".$project."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "#..............................................................................\n");
	fwrite($logOutput, "Running 'php/project.paired_WGseq.install_2.php'.\n");
	fwrite($logOutput, "Variables passed via command-line from 'php/project.paired_WGseq.install_1.php' :\n");
	fwrite($logOutput, "\tfileName         = '".$fileName."'\n");
	fwrite($logOutput, "\tuser             = '".$user."'\n");
	fwrite($logOutput, "\tproject          = '".$project."'\n");
	fwrite($logOutput, "#============================================================================== 3\n");

// Manage condensed log file.
	$condensedLogOutputName = "../users/".$user."/projects/".$project."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'a');
//	fclose($condensedLogOutput);

// Generate 'datafiles.txt' file containing: name of all data files.
// Identify format of uploaded file and decompress as needed (*.ZIP; *.GZ).
	$outputName = "../users/".$user."/projects/".$project."/datafiles.txt";
	$output     = fopen($outputName, 'w');
	$fileNames = explode(",", $fileName);
	fwrite($logOutput, "\tGenerate 'datafiles.txt' and decompress uploaded archives.\n");
	foreach ($fileNames as $key=>$name) {
		$projectPath = "../users/".$user."/projects/".$project."/";
		$name        = str_replace("\\", ",", $name);
		$ext         = strtolower(pathinfo($name, PATHINFO_EXTENSION));
		$filename    = strtolower(pathinfo($name, PATHINFO_FILENAME));
		fwrite($logOutput, "\tFile ".$key."\n");
		fwrite($logOutput, "\t\tDatafile   : '$name'.\n");
		fwrite($logOutput, "\t\tFilename   : '$filename.'.\n");
		fwrite($logOutput, "\t\tExtension  : '$ext'.\n");
		fwrite($logOutput, "\t\tPath       : '$projectPath'.\n");
		// Process the uploaded file.
		$paired = process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput);
		// formatting.
		if ($key < count(fileNames)-1) {
			fwrite($output,"\n");
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
	$system_call_string = "sh ../sh/project.paired_WGseq.install_3.sh ".$user." ".$project." > /dev/null &";
	system($system_call_string);
	fclose($condensedLogOutput);
	fclose($logOutput);
?>
