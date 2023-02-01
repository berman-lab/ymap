<?php
	session_start();
        error_reporting(E_ALL);
        require_once 'constants.php';
	require_once 'POST_validation.php';
        ini_set('display_errors', 1);

        // If the user is not logged on, redirect to login page.
        if(!isset($_SESSION['logged_on'])){
		session_destroy();
                header('Location: .');
        }

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Sanitize input strings.
	$genome = sanitize_POST("genome");
	$key    = sanitize_POST("key");
	$status = sanitize_POST("status");

	// Confirm if requested genome exists.
	$genome_dir = "users/".$user."/genomes/".$genome;
	if (!is_dir($genome_dir)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
	}

?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">
<head>
<style type="text/css">
	body {font-family: arial;}
	.tab {
		margin-left:    1cm;
	}
	.clock {
		float:          left;
		margin-right:   0.25cm;
	}
</style>
</head>
<?php
	// increment clock animation...
	$status   = ($status + 1) % 12;
	if ($status == 0) {          $clock = "<img src=\"images/12.png\" alt-text=\"12\" class=\"clock\" >";
	} else if ($status == 1) {   $clock = "<img src=\"images/01.png\" alt-text=\"01\" class=\"clock\" >";
	} else if ($status == 2) {   $clock = "<img src=\"images/02.png\" alt-text=\"02\" class=\"clock\" >";
	} else if ($status == 3) {   $clock = "<img src=\"images/03.png\" alt-text=\"03\" class=\"clock\" >";
	} else if ($status == 4) {   $clock = "<img src=\"images/04.png\" alt-text=\"04\" class=\"clock\" >";
	} else if ($status == 5) {   $clock = "<img src=\"images/05.png\" alt-text=\"05\" class=\"clock\" >";
	} else if ($status == 6) {   $clock = "<img src=\"images/06.png\" alt-text=\"06\" class=\"clock\" >";
	} else if ($status == 7) {   $clock = "<img src=\"images/07.png\" alt-text=\"07\" class=\"clock\" >";
	} else if ($status == 8) {   $clock = "<img src=\"images/08.png\" alt-text=\"08\" class=\"clock\" >";
	} else if ($status == 9) {   $clock = "<img src=\"images/09.png\" alt-text=\"09\" class=\"clock\" >";
	} else if ($status == 10) {  $clock = "<img src=\"images/10.png\" alt-text=\"10\" class=\"clock\" >";
	} else if ($status == 11) {  $clock = "<img src=\"images/11.png\" alt-text=\"11\" class=\"clock\" >";
	} else {                     $clock = "[ * ]";
	}

	$sizeFile_1   = $genome_dir."/upload_size_1.txt";
	$handle       = fopen($sizeFile_1,'r');
	$sizeString_1 = trim(fgets($handle));
	fclose($handle);
	if ($sizeString_1 !== "") {
		echo "\n<script type='text/javascript'>\n";
		echo "parent.parent.update_genome_file_size('".$key."','".$sizeString_1."');";
		echo "\n</script>\n";
	}

	if (file_exists($genome_dir."/complete.txt")) {
		?>
		<html>
		<body onload = "parent.parent.update_genome_label_color('<?php echo $key; ?>','#00AA00'); parent.parent.update_genome_remove_iframe('<?php echo $key; ?>');">
		</body>
		</html>
		<?php
	} else if (file_exists($genome_dir."/working.txt")) {
		// Load start time from 'working.txt'
		$startTimeStamp = file_get_contents($genome_dir."/working.txt");
		$startTime      = strtotime($startTimeStamp);
		$currentTime    = time();
		$intervalTime   = $currentTime - $startTime;
		if ($intervalTime > 60*60*24) { // likely error.
			?>
			<BODY onload = "parent.parent.resize_genome('<?php echo $key; ?>', 100);" class="tab">
				<font color="red" size="2"><b>[Error]</b></font><?php echo " &nbsp; &nbsp; ".$clock; ?><br>
				<font size="2">Processing of genome data has taken longer than expected and might be stalled.<br>Contact the admin through the "System" tab with details and they will check on the job.</font>
			</BODY>
			</HTML>
			<?php
		} else {
			// Load last line from "condensed_log.txt" file.
			$condensedLog      = explode("\n", trim(file_get_contents($genome_dir."/condensed_log.txt")));
			$condensedLogEntry = $condensedLog[count($condensedLog)-1];
			?>
			<script type="text/javascript">
			var user    = "<?php echo $user; ?>";
			var genome  = "<?php echo $genome; ?>";
			var key     = "<?php echo $key; ?>";
			var status  = "<?php echo $status; ?>";
			reload_page=function() {
			<?php
				// Make a form to generate a form to POST information to pass along to page reloads, auto-triggered by form submit.
				echo "\tvar autoSubmitForm = document.createElement('form');\n";
				echo "\t\tautoSubmitForm.setAttribute('method','post');\n";
				echo "\t\tautoSubmitForm.setAttribute('action','genome.working_server.php');\n";
				echo "\tvar input2 = document.createElement('input');\n";
				echo "\t\tinput2.setAttribute('type','hidden');\n";
				echo "\t\tinput2.setAttribute('name','key');\n";
				echo "\t\tinput2.setAttribute('value',key);\n";
				echo "\t\tautoSubmitForm.appendChild(input2);\n";
				echo "\tvar input2 = document.createElement('input');\n";
				echo "\t\tinput2.setAttribute('type','hidden');\n";
				echo "\t\tinput2.setAttribute('name','user');\n";
				echo "\t\tinput2.setAttribute('value',user);\n";
				echo "\t\tautoSubmitForm.appendChild(input2);\n";
				echo "\tvar input3 = document.createElement('input');\n";
				echo "\t\tinput3.setAttribute('type','hidden');\n";
				echo "\t\tinput3.setAttribute('name','genome');\n";
				echo "\t\tinput3.setAttribute('value',genome);\n";
				echo "\t\tautoSubmitForm.appendChild(input3);\n";
				echo "\tvar input4 = document.createElement('input');\n";
				echo "\t\tinput4.setAttribute('type','hidden');\n";
				echo "\t\tinput4.setAttribute('name','status');\n";
				echo "\t\tinput4.setAttribute('value',status);\n";
				echo "\t\tautoSubmitForm.appendChild(input4);\n";
				echo "\t\tdocument.body.appendChild(autoSubmitForm);\n";
				echo "\tautoSubmitForm.submit();\n";
				?>
			}
			// Initiate recurrent call to reload_page function, which depends upon genome status.
			var internalIntervalID = window.setInterval(reload_page, 3000);
			</script>
			<BODY onload = "parent.parent.resize_genome('<?PHP echo $key; ?>', 38);" class="tab">
				<font color="red" size="2"><b>[Processing uploaded data.]</b></font><?php echo " &nbsp; &nbsp; ".$clock; ?><br>
				<font size="2"><?php
				echo $condensedLogEntry;
				?></font>
			</BODY>
			</HTML>
			<?php
		}
	}
?>
