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
//	$fileName = $argv[1];
//	$user     = $argv[2];
//	$project  = $argv[3];
//	$key      = $argv[4];

	$fileName = "";
	$user     = "darren";
	$project  = "morlurie_30-26";
	$key      = "0";

// Initialize log file.
	$logOutputName = $directory."users/".$user."/projects/".$project."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "#..............................................................................\n");
	fwrite($logOutput, "Running 'php/project.single_WGseq.install_2.php'.\n");
	fwrite($logOutput, "Variables passed via command-line from 'php/project.single_WGseq.install_1.php' :\n");
	fwrite($logOutput, "\tfileName = '".$fileName."'\n");
	fwrite($logOutput, "\tuser     = '".$user."'\n");
	fwrite($logOutput, "\tproject  = '".$project."'\n");
	fwrite($logOutput, "\tkey      = '".$key."'\n");
	fwrite($logOutput, "#============================================================================== 3\n");

// Manage condensed log file.
	$condensedLogOutputName = $directory."users/".$user."/projects/".$project."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'a');
//	fclose($condensedLogOutput);

	// Final install functions are in shell script.
	fwrite($logOutput, "Passing control to : 'sh/project.single_WGseq.install_3.sh'\n");
 	$system_call_string = "sh ../sh/project.single_WGseq.install_3.sh ".$user." ".$project." ".$directory." > /dev/null &";
	system($system_call_string);

	fclose($condensedLogOutput);
?>
