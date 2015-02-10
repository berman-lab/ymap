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
		.tab {
			margin:        40px;
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
	$user            = filter_input(INPUT_POST, "user",            FILTER_SANITIZE_STRING);
	$project         = filter_input(INPUT_POST, "project",         FILTER_SANITIZE_STRING);

// Initialize 'process_log.txt' file.
	$outputLogName = "../users/".$user."/projects/".$project."/process_log.txt";
	$outputLog     = fopen($outputLogName, 'w');
    fwrite($outputLog, "Process_log.txt initialized.\n");

// Call Matlab to process SnpCgh data file into final figure.
	$designDefinition   = "design1";
	$inputFile          = "../../users/".$user."/projects/".$project."/".$fileName;
	$headerRows         = 46;
	$colNames           = 1;
	$colCh1             = 4;
	$colCh2             = 5;
	$colRatio           = 6;
	$colLog2ratio       = 7;
	$phasingData        = "cal_paper";
	$ploidyFileContents = file_get_contents("../users/".$user."/projects/".$project."/ploidy.txt");
		$ploidyStrings  = explode("\n", $ploidyFileContents);
		$ploidyEstimate = $ploidyStrings[0];
		$ploidyBase     = $ploidyStrings[1];

	$imageFormat        = "png";
	$projectName        = $project;
	$workingDir         = "../../users/".$user."/projects/".$project."/";
	$show_MRS           = file_get_contents("../users/".$user."/projects/".$project."/showAnnotations.txt");

	$outputName         = "../users/".$user."/projects/".$project."/processing.m";
	$dirBase            = "../users/".$user."/projects/".$project."/";
	$output             = fopen($outputName, 'w');
	$outputString       =  "function [] = processing()\n";
	$outputString      .= "\tdiary('matlab.process_log.txt');\n";
	$outputString      .= "\tcd ../../../../Matlab/SnpCgh_array;\n";

	$outputString      .= "\tprocess_main_reinstall('".$phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";
	$outputString      .= "\tprocess_main_CNV_only('".$designDefinition."','".$inputFile."','".$headerRows."','".$colNames."','".$colCh1."','".$colCh2."','".$colRatio."','";
	$outputString      .=    $colLog2ratio."','".$phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";
	$outputString      .= "\tprocess_main_SNP_only('".$designDefinition."','".$inputFile."','".$headerRows."','".$colNames."','".$colCh1."','".$colCh2."','".$colRatio."','";
	$outputString      .=    $colLog2ratio."','".$phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";

	$outputString      .= "\t\n";
	$outputString      .= "end";
	fwrite($output, $outputString);
	fclose($output);
	chmod($outputName,0644);
	fwrite($outputLog, "Matlab running script 'processing.m' generated.\n");

	fwrite($outputLog, "Current Path = '".getcwd()."'\n");
	fwrite($outputLog, "Calling Matlab :\n");
	$system_call_string_2 = "matlab -nosplash -nodesktop 'run ../users/".$user."/projects/".$project."/processing.m' > /dev/null &";
	fwrite($outputLog, "\t\"".$system_call_string_2."\"\n\n");
	fclose($outputLog);
	system($system_call_string_2);
?>
</script>
</BODY>
</HTML>
