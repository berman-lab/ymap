<?php
	session_start();
	error_reporting(E_ALL);
	require_once '../constants.php';
	require_once '../POST_validation.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: ../');
	}

	// pull strings from session.
	$user     = $_SESSION['user'];
	$fileName = $_SESSION['fileName'];
	$project  = $_SESSION['project'];
	$key      = $_SESSION['key'];

	$project_dir = "../users/".$user."/projects/".$project;
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
// Initialize 'process_log.txt' file.
	$outputLogName   = $project_dir."/process_log.txt";
	$outputLog       = fopen($outputLogName, 'w');
	fwrite($outputLog, "Process_log.txt initialized.\n");

// Generate 'working.txt' file to let pipeline know processing is started.
	$outputName      = $project_dir."/working.txt";
	$output          = fopen($outputName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($output, $startTimeString);
	fclose($output);
	chmod($outputName,0644);
	fwrite($outputLog, "'working.txt' file generated.\n");

// Generate 'upload_size.txt' file to contain the size of the uploaded file (irrespective of format) for display in "Manage Datasets" tab.
	$outputName      = $project_dir."/upload_size_1.txt";
	$output          = fopen($outputName, 'w');
	$fileSizeString  = filesize($project_dir."/".$fileName);
	fwrite($output, $fileSizeString);
	fclose($output);
	chmod($outputName,0644);
	fwrite($outputLog, "\tGenerated 'upload_size_1.txt' file.\n");

// Generate 'datafile1.txt' file containing: name of first data file.
	$outputName = $project_dir."/datafile1.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, $fileName);
	fclose($output);
	chmod($outputName,0644);
	fwrite($outputLog, "'datafile1.txt' file generated.\n");

// Call Matlab to process SnpCgh data file into final figure.
	$designDefinition   = "design1";
	$inputFile          = $project_dir."/".$fileName;
	$headerRows         = 46;
	$colNames           = 1;
	$colCh1             = 4;
	$colCh2             = 5;
	$colRatio           = 6;
	$colLog2ratio       = 7;
	$phasingData        = "cal_paper";
	$ploidyFileContents = file_get_contents($project_dir."/ploidy.txt");
		$ploidyStrings  = explode("\n", $ploidyFileContents);
		$ploidyEstimate = $ploidyStrings[0];
		$ploidyBase     = $ploidyStrings[1];

	$imageFormat        = "png";
	$projectName        = $project;
	$workingDir         = $project_dir."/";
	$show_MRS           = file_get_contents($project_dir."/showAnnotations.txt");

	$outputName         = $project_dir."/processing.m";
	$dirBase            = $project_dir."/";
	$output             = fopen($outputName, 'w');
	$outputString       =  "function [] = processing()\n";
	$outputString      .= "\tdiary('matlab.process_log.txt');\n";
	$outputString      .= "\tcd ../../../../scripts_SnpCghArray;\n";

	// Log status to process_log.txt file in project directory.
	$outputString      .= "\tnew_fid = fopen(".$project_dir."/process_log.txt','a');\n";
	$outputString      .= "\tfprintf(new_fid,'Starting \"process_main.m\".\\n');\n";

	$outputString      .= "\tprocess_main('".$designDefinition."','".$inputFile."','".$headerRows."','".$colNames."','".$colCh1."','".$colCh2."','".$colRatio."','".$colLog2ratio."','";
	$outputString      .=    $phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";

	// Log status to process_log.txt file in project directory.
	$outputString      .= "\tfprintf(new_fid,'Starting \"process_main_CNV_only.m\".\\n');\n";

	$outputString      .= "\tprocess_main_CNV_only('".$designDefinition."','".$inputFile."','".$headerRows."','".$colNames."','".$colCh1."','".$colCh2."','".$colRatio."','";
	$outputString      .=    $colLog2ratio."','".$phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";

	// Log status to process_log.txt file in project directory.
	$outputString      .= "\tfprintf(new_fid,'Starting \"process_main_SNP_only.m\".\\n');\n";

	$outputString      .= "\tprocess_main_SNP_only('".$designDefinition."','".$inputFile."','".$headerRows."','".$colNames."','".$colCh1."','".$colCh2."','".$colRatio."','";
	$outputString      .=    $colLog2ratio."','".$phasingData."','".$ploidyEstimate."','".$ploidyBase."','".$imageFormat."','".$projectName."','".$workingDir."','".$show_MRS."');\n";

	// Log status to process_log.txt file in project directory.
	$outputString      .= "\tfclose(new_fid);\n";

	$outputString      .= "end";
	fwrite($output, $outputString);
	fclose($output);
	chmod($outputName,0644);
	fwrite($outputLog, "Matlab running script 'processing.m' generated.\n");

	// Running pre-processing for Matlab script, because running it within Matlab doesn't work well:
	// http://stackoverflow.com/questions/29451735/matlab-system-function-not-running-the-command-and-not-returning-any-values

	// Standardize end-of-line characters to '\n':
	system('perl -pi -e "s/\r\n/\n/g" ' . $inputFile); // \r\n => \n
	system('perl -pi -e "s/\r/\n/g" ' . $inputFile); // \r   => \n

	/*
	Make CGH/SNP/MLST probe data subfiles.

	CGH probes  : "CGHv1_Ca_..."
	SNP probes  : "SNPv1_Ca_..."
	MLST probes : "MLSTv1_..."
	*/
	$cghRowsFile = $workingDir . 'CGH_rows.xls';
	$snpRowsFile = $workingDir . 'SNP_rows.xls';
	$mlstRowsFile = $workingDir . 'MLST_rows.xls';
	system('grep "CGHv1_Ca_\|CGH_Ca_" ' . $inputFile . ' > ' . $cghRowsFile);
	system('grep "SNPv1_Ca_\|SNP_Ca_" ' . $inputFile . ' > ' . $snpRowsFile);
	system('grep "MLSTv1_" ' . $inputFile . ' > ' . $mlstRowsFile);

	/*
	Sort subfiles by probe name

	-b    : ignore leading spaces in column entries.
	-k1,1 : sort on columns from 1 to 1.
	-o    : output file.
	*/
	foreach (array($cghRowsFile, $snpRowsFile, $mlstRowsFile) as $fileToProcess)
	{
		chmod($fileToProcess, 0644);
		$tmpFile = "$fileToProcess.tmp";
		system("sort -b -k1,1 $fileToProcess -o $tmpFile");
		system("mv $tmpFile $fileToProcess");
	}

	fwrite($outputLog, "Current Path = '".getcwd()."'\n");
	fwrite($outputLog, "Calling Matlab :\n");
	$system_call_string_2 = ". ../local_installed_programs.sh && $matlab_exec -nosplash -nodesktop -r \'run ".$project_dir."/processing.m\' > /dev/null &";
	fwrite($outputLog, "\t\"".$system_call_string_2."\"\n\n");
	fclose($outputLog);
	system($system_call_string_2);
?>
<BODY>
<font size="2" color="red">Upload complete; processing...</font><br>
<script type="text/javascript">
	var user    = "<?php echo $user;    ?>";
	var project = "<?php echo $project; ?>";
	var key     = "<?php echo $key;     ?>";
	var status  = "0";
	// construct and submit form to move on to "project.working_server.php";
	var autoSubmitForm = document.createElement("form");
		autoSubmitForm.setAttribute("method","post");
		autoSubmitForm.setAttribute("action","../project.working_server.php");
	var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","key");
		input2.setAttribute("value",key);
		autoSubmitForm.appendChild(input2);
	var input2 = document.createElement("input");
		input2.setAttribute("type","hidden");
		input2.setAttribute("name","user");
		input2.setAttribute("value",user);
		autoSubmitForm.appendChild(input2);
	var input3 = document.createElement("input");
		input3.setAttribute("type","hidden");
		input3.setAttribute("name","project");
		input3.setAttribute("value",project);
		autoSubmitForm.appendChild(input3);
	var input4 = document.createElement("input");
		input4.setAttribute("type","hidden");
		input4.setAttribute("name","status");
		input4.setAttribute("value",status);
	document.body.appendChild(autoSubmitForm);
	autoSubmitForm.submit();
</script>
</BODY>
</HTML>
