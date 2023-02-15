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
	require_once '../../constants.php';
	include_once '../../process_input_files.php';

// Deal with passed variables.
	$fileName = $argv[1];
	$user     = $argv[2];
	$project  = $argv[3];

// Initialize log file.
	$logOutputName = "../../users/".$user."/projects/".$project."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "#..............................................................................\n");
	fwrite($logOutput, "Running 'scripts_seqModules/scripts_ddRADseq/project.single_ddRADseq.install_2.php'.\n");
	fwrite($logOutput, "Variables passed via command-line from 'scripts_seqModules/scripts_ddRADseq/project.single_ddRADseq.install_1.php' :\n");
	fwrite($logOutput, "\tfileName = '".$fileName."'\n");
	fwrite($logOutput, "\tuser     = '".$user."'\n");
	fwrite($logOutput, "\tproject  = '".$project."'\n");
	fwrite($logOutput, "#============================================================================== 3\n");

// Manage condensed log file.
	$condensedLogOutputName = "../../users/".$user."/projects/".$project."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'a');
//	fclose($condensedLogOutput);

// Generate 'datafiles.txt' file containing: name of all data files.
// Identify format of uploaded file and decompress as needed (*.ZIP; *.GZ).
	$outputName = "../../users/".$user."/projects/".$project."/datafiles.txt";
	$output     = fopen($outputName, 'w');
	$fileNames  = explode(",", $fileName);
	fwrite($logOutput, "\tGenerate 'datafiles.txt' and decompress uploaded archives.\n");
	$paired     = 0;
	foreach ($fileNames as $key=>$name) {
		$projectPath = "../../users/".$user."/projects/".$project."/";
		$name        = str_replace("\\", ",", $name);
		//rename($projectPath.$name,$projectPath.strtolower($name));
		//$name        = strtolower($name);
		$ext         = strtolower(pathinfo($name, PATHINFO_EXTENSION));
		$filename    = strtolower(pathinfo($name, PATHINFO_FILENAME));
		fwrite($logOutput, "\tFile ".$key."\n");
		fwrite($logOutput, "\t\tDatafile  : '$name'.\n");
		fwrite($logOutput, "\t\tFilename  : '$filename.'.\n");
		fwrite($logOutput, "\t\tExtension : '$ext'.\n");
		fwrite($logOutput, "\t\tPath      : '$projectPath'.\n");

		// Generate 'upload_size.txt' file to contain the size of the uploaded file (irrespective of format) for display in "Manage Datasets" tab.
        $output2Name    = $projectPath."upload_size_1.txt";
        $output2        = fopen($output2Name, 'w');
        $fileSizeString = filesize($projectPath.$name);
        fwrite($output2, $fileSizeString);
        fclose($output2);
        chmod($output2Name,0755);
        fwrite($logOutput, "\tGenerated 'upload_size_1.txt' file.\n");

		// Process the uploaded file.
		$paired = process_input_files($ext,$name,$projectPath,$key,$user,$project,$output, $condensedLogOutput,$logOutput);
		// formatting.
		if ($key < count(fileNames)-1) {
			fwrite($output,"\n");
		}
	}
	fclose($output);
	chmod($outputName,0755);
	fwrite($logOutput, "Completed 'datafiles.txt' file.\n");

	// Final install functions are in shell script.
	if ($paired == 1) {
		fwrite($logOutput, "Passing control to : 'scripts_seqModules/scripts_ddRADseq/project.paired_ddRADseq.install_3.sh'\n");
		fwrite($logOutput, "\t\tPaired-end reads being processed.\n");
		fwrite($logOutput, "\t\tCurrent directory = '".getcwd()."'\n" );
		$system_call_string = "sh project.paired_ddRADseq.install_3.sh ".$user." ".$project." > /dev/null &";
	} else {
		fwrite($logOutput, "Passing control to : 'scripts_seqModules/scripts_ddRADseq/project.single_ddRADseq.install_3.sh'\n");
		fwrite($logOutput, "\t\tSingle-end reads being processed.\n");
		fwrite($logOutput, "\t\tCurrent directory = '".getcwd()."'\n" );
	 	$system_call_string = "sh project.single_ddRADseq.install_3.sh ".$user." ".$project." > /dev/null &";
	}
	fwrite($logOutput, "\t\tSystem call string = '".$system_call_string."'\n");

	system($system_call_string);
	fclose($condensedLogOutput);
	fclose($logOutput);
?>
