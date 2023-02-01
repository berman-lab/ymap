<?php
	// Attempt to setup php for disconnecting from web browser.
	ob_end_clean();
	header("Connection: close");
	ob_start();

	session_start();
	error_reporting(E_ALL);
        require_once '../constants.php';
        ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
		session_destroy();
                header('Location: user.login.php');
        }

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Sanitize input strings.
	$key    = trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING));
	$key    = str_replace(" ","_",$key);
	$key    = preg_replace("/[\s\W]+/", "", $key);

	// Load genome string from session.
	$genome   = $_SESSION['genome_'.$key];
	$genome_dir = "../users/".$user."/genomes/".$genome;

	// Sanitize input strings.
	$fileName = trim(filter_input(INPUT_POST, "fileName", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$fileName = str_replace(" ","_",$fileName);					// convert any spaces to underlines.
	$fileName = preg_replace("/[\s\W]+/", "", $fileName);				// remove everything but alphanumeric characters and underlines.

	// load PHP function to process input files.
	include_once 'process_input_files.genome.php';
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
<title>Install genome into pipeline.</title>
</HEAD>
<?php
	// Open 'process_log.txt' file.
	$logOutputName = $genome_dir."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'scripts_genomes/genome.install_5.php'.\n");
	fwrite($logOutput, "\tkey     :'".$key."'\n");
	if ($filename == '') {
		fwrite($logOutput, "\tfileName: [No chromosome features file loaded].\n");
	} else {
		fwrite($logOutput, "\tfileName:'".$fileName."'\n");
	}

	// Generate 'working.txt' to tell main page that genome installation is in process.
	fwrite($logOutput, "\tGenerating 'working.txt' file.\n");
	$outputName      = $genome_dir."/working.txt";
	$output          = fopen($outputName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($output, $startTimeString);
	fclose($output);

	if ($fileName == '') {
		// No uploaded chromosome features file.
		fwrite($logOutput, "\tNo chromosome features file uploaded.\n");
	} else {
		// Process uploaded file.
		$genomePath  = $genome_dir."/";
		$name        = str_replace("\\", ",", $fileName);
		rename($genomePath.$name,$genomePath.strtolower($name));
		$name        = strtolower($name);
		$ext         = strtolower(pathinfo($name, PATHINFO_EXTENSION));
		$filename    = strtolower(pathinfo($name, PATHINFO_FILENAME));
		$newName    = "chromosome_features.txt";

		// Generate 'upload_size.txt' file to contain the size of the uploaded file (irrespective of format) for display in "Manage Datasets" tab.
		$output2Name    = $genomePath."upload_size_2.txt";
		$output2        = fopen($output2Name, 'w');
		$fileSizeString = filesize($genomePath.$name);
		fwrite($output2, $fileSizeString);
		fclose($output2);
		chmod($output2Name,0755);
		fwrite($logOutput, "\tGenerated 'upload_size_1.txt' file.\n");

		// Process the uploaded file.
		process_input_files_genome($ext,$name,$genomePath,$key,$user,$genome,$output, $condensedLogOutput,$logOutput, $newName);
		$fileName = $newName;
		fwrite($logOutput, "\tProcessed chromosome features file.\n");
	}

	// Final install functions are in shell script.
	$system_call_string = "sh genome.install_6.sh ".$user." ".$genome." > /dev/null &";
	system($system_call_string);
	fwrite($logOutput, "\tCurrent Directory  = '".getcwd()."'\n");
	fwrite($logOutput, "\tSystem Call String = '".$system_call_string."'\n");

	// Debugging output of all variables.
	// print_r($GLOBALS);
	// exit;

// The following section is to trigger the interface componant which shows status of the process, while the process has already been spawned off above.
?>
<BODY onload = "parent.parent.resize_genome('<?php echo $key; ?>', 40);">
	<font color="red">[Installation in process.]</font>
	<font size="2" color="red">Upload complete; processing...</font><br>
	<script type="text/javascript">
		var user   = "<?php echo $user;   ?>";
		var genome = "<?php echo $genome; ?>";
		var key    = "<?php echo $key;    ?>";
		var status = "0";
		// construct and submit form to move on to "project.working_server.2.php";
		var autoSubmitForm = document.createElement("form");
			autoSubmitForm.setAttribute("method","post");
			autoSubmitForm.setAttribute("action","../genome.working_server.php");
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
			input3.setAttribute("name","genome");
			input3.setAttribute("value",genome);
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
<?php
// process_log.txt output.
	fwrite($logOutput, "\t'scripts_genomes/genome.install_5.php' has completed.\n");
	fclose($logOutput);

// Initialize 'condensed_log.txt' file.
	$condensedLogOutputName = $genome_dir."/condensed_log.txt";
	$condensedLogOutput     = fopen($condensedLogOutputName, 'w');
	fwrite($condensedLogOutput, "Process started.\n");
	fclose($condensedLogOutput);
	chmod($outputName,0755);
?>
