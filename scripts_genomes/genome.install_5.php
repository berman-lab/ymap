<?php
	// Attempt to setup php for disconnecting from web browser.
	ob_end_clean();
	header("Connection: close");
	ob_start();

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
<title>Install genome into pipeline.</title>
</HEAD>
<?php
    require_once '../constants.php';

	$user             = $_SESSION['user'];
	$key              = filter_input(INPUT_POST, "key",      FILTER_SANITIZE_STRING);
	$fileName         = filter_input(INPUT_POST, "fileName", FILTER_SANITIZE_STRING);
	$genome           = $_SESSION['genome_'.$key];

//	$chr_count        = $_SESSION['chr_count_'.$key];
//	$chr_used_count   = $_SESSION['chr_used_count_'.$key];
//	$chr_lengths      = $_SESSION['chr_lengths_'.$key];
//	$chr_names        = $_SESSION['chr_names_'.$key];
//	$chr_draws        = $_SESSION['chr_draws_'.$key];
//	$chr_shortNames   = $_SESSION['chr_shortNames_'.$key];
//	$chr_cenStarts    = $_SESSION['chr_cenStarts_'.$key];
//	$chr_cenEnds      = $_SESSION['chr_cenEnds_'.$key];
//	$rDNA_chr         = $_SESSION['rDNA_chr_'.$key];
//	$rDNA_start       = $_SESSION['rDNA_start_'.$key];
//	$rDNA_end         = $_SESSION['rDNA_end_'.$key];
//	$ploidyDefault    = $_SESSION['ploidyDefault_'.$key];
//	$annotation_count = $_SESSION['annotation_count_'.$key];

// Open 'process_log.txt' file.
	$logOutputName = "../users/".$user."/genomes/".$genome."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'scripts_genomes/genome.install_5.php'.\n");
	fwrite($logOutput, "\tkey     :'".$key."'\n");
	fwrite($logOutput, "\tfileName:'".$fileName."'\n");

// Generate 'working.txt' to tell main page that genome installation is in process.
	fwrite($logOutput, "\tGenerating 'working.txt' file.\n");
	$outputName      = "../users/".$user."/genomes/".$genome."/working.txt";
	$output          = fopen($outputName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($output, $startTimeString);
	fclose($output);

// process_log.txt output.
	fwrite($logOutput, "\t'scripts_genomes/genome.install_5.php' has completed.\n");

// Standardize name of uploaded file.
	if ($fileName == '') {
	} else {
		$oldName    = $fileName;
		$newName    = "chromosome_features.txt";
		$currentDir = getcwd();
		$genomePath = "../users/".$user."/genomes/".$genome."/";
		chdir($genomePath);
		rename($oldName, $newName);
		chdir($currentDir);
	}

// Final install functions are in shell script.
	$system_call_string = "sh ../sh/genome.install_6.sh ".$user." ".$genome." > /dev/null &";
	system($system_call_string);
	fwrite($logOutput, getcwd()."\n");
	fwrite($logOutput, "'".$system_call_string."'\n");

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
		autoSubmitForm.submit();
	</script>
</BODY>
</HTML>
