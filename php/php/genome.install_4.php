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
    require_once 'constants.php';

	$user             = $_SESSION['user'];
	$key              = filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING);
	$genome           = $_SESSION['genome_'.$key];

// Open 'process_log.txt' file.
	$logOutputName = $directory."users/".$user."/genomes/".$genome."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'php/genome.install_4.php'.\n");

// process POST data.
	fwrite($logOutput, "\tProcessing POST data containing genome specific information.\n");
	$headerLineCount = filter_input(INPUT_POST, "headerLineCount", FILTER_SANITIZE_STRING);
	$col_chrID       = filter_input(INPUT_POST, "col_chrID",       FILTER_SANITIZE_STRING);
	$col_startBP     = filter_input(INPUT_POST, "col_startBP",     FILTER_SANITIZE_STRING);
	$col_endBP       = filter_input(INPUT_POST, "col_endBP",       FILTER_SANITIZE_STRING);

// Generate 'expression.txt' to tell main page that genome installation is in process.
	fwrite($logOutput, "\tGenerating 'expresison.txt' file.\n");
	$outputName = $directory."users/".$user."/genomes/".$genome."/expression.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, $headerLineCount."\n");
	fwrite($output, $col_chrID."\n");
	fwrite($output, $col_startBP."\n");
	fwrite($output, $col_endBP."\n");
	fclose($output);

// process_log.txt output.
	fwrite($logOutput, "\t'php/genome.install_4.php' has passed processing off to 'sh/genome.install_5.sh'.\n");

	// Final install functions are in shell script.
	$system_call_string = "sh ../sh/genome.install_5.sh ".$user." ".$genome." ".$directory." > /dev/null &";
	system($system_call_string);
?>
<BODY onload = "parent.resize_iframe('<?php echo $key; ?>', 40);">
	<div id='frameContainer'></div>
	<script type="text/javascript">
		<?php
		echo "var el_g = document.getElementById('frameContainer');\n\t\t";
		echo "el_g.innerHTML='<iframe id=\"g_new\" name=\"g_new\" class=\"upload\" style=\"height:38px\" src=\"../uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
		echo "frames['g_new'].display_string = new Array();\n\t\t";
		echo "frames['g_new'].display_string[0] = \"Add : chromosome feature file...\";\n\t\t";
		echo "frames['g_new'].target_dir        = \"".$directory."users/".$user."/genomes/".$genome."/\";\n\t\t";
		echo "frames['g_new'].conclusion_script = \"".$url."php/genome.install_5.php\";\n\t\t";
		echo "frames['g_new'].user              = \"".$user."\";\n\t\t";
		echo "frames['g_new'].genome            = \"".$genome."\";\n\t\t";
		echo "frames['g_new'].key               = \"".$key."\";\n\t";
?></script>
</BODY>
</HTML>
