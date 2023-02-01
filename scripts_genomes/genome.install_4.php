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
	$genome = trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$genome = str_replace(" ","_",$genome);						// convert any spaces to underlines.
	$genome = preg_replace("/[\s\W]+/", "", $genome);				// remove everything but alphanumeric characters and underlines.
	$key    = trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING));
	$key    = str_replace(" ","_",$key);
	$key    = preg_replace("/[\s\W]+/", "", $key);

	// Confirm if requested genome exists.
	$genome_dir = "../users/".$user."/genomes/".$genome;
	if (!is_dir($genome_dir)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}
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
	fwrite($logOutput, "Running 'scripts_genomes/genome.install_4.php'.\n");

// process POST data.
	fwrite($logOutput, "\tProcessing POST data containing genome specific information.\n");
	$headerLineCount = filter_input(INPUT_POST, "headerLineCount", FILTER_SANITIZE_STRING);
	$col_chrID       = filter_input(INPUT_POST, "col_chrID",       FILTER_SANITIZE_STRING);
	$col_startBP     = filter_input(INPUT_POST, "col_startBP",     FILTER_SANITIZE_STRING);
	$col_endBP       = filter_input(INPUT_POST, "col_endBP",       FILTER_SANITIZE_STRING);

// Generate 'expression.txt' to tell main page that genome expression analysis installation is in process.
	fwrite($logOutput, "\tGenerating 'expresison.txt' file.\n");
	$outputName = $genome_dir."/expression.txt";
	$output     = fopen($outputName, 'w');
	fwrite($output, $headerLineCount."\n");
	fwrite($output, $col_chrID."\n");
	fwrite($output, $col_startBP."\n");
	fwrite($output, $col_endBP."\n");
	fclose($output);

// process_log.txt output.
	fwrite($logOutput, "\t'scripts_genomes/genome.install_4.php' has completed'.\n");

?>
<BODY onload = "parent.parent.resize_genome('<?php echo $key; ?>', 40);">
	<div id='frameContainer'></div>
	<script type="text/javascript">
		<?php
		echo "var el_g = document.getElementById('frameContainer');\n\t\t";
		echo "el_g.innerHTML='<iframe id=\"g_new\" name=\"g_new\" class=\"upload\" style=\"height:38px\" src=\"../uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
		echo "frames['g_new'].display_string    = new Array();\n\t\t";
		echo "frames['g_new'].display_string[0] = \"Add : chromosome feature TAB file...\";\n\t\t";
		echo "frames['g_new'].target_dir        = \"../".$genome_dir."/\";\n\t\t";
		echo "frames['g_new'].conclusion_script = \"scripts_genomes/genome.install_5.php\";\n\t\t";
		echo "frames['g_new'].user              = \"".$user."\";\n\t\t";
		echo "frames['g_new'].genome            = \"".$genome."\";\n\t\t";
		echo "frames['g_new'].key               = \"".$key."\";\n\t";
?></script>
</BODY>
</HTML>
