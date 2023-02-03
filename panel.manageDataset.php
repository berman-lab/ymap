<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>
	Install SNP/CGH microarray and next generation sequence datasets to, or delete datasets from, your user account.
	<font size='2'>
        <p>
            <b>IMPORTANT:</b>
            <ol>
                <li>Filenames should only have alphanumeric characters (letters, numbers, underscores and dashes) in their names (no spaces or other special characters!).</li>
                <li>.fq and .fq.gz file extensions are not supported - rename them to .fastq or .fastq.gz (respectively).</li>
                <li>Use Chrome for uploads (this will allow resuming them if they are interrupted).</li>
            </ol>
        </p>
		<p>
			Is your upload stuck? To resume it, do the following:
			<ol>
				<li>When all other uploads are finished, refresh the page.</li>
				<li>Add the exact same files as you have previously to the project that is stuck.</li>
				<li>Press the "upload" button.</li>
			</ol>
		</p>
	</font>
</font>
<?php
	if (isset($_SESSION['logged_on']))
	{
		require_once 'sharedFunctions.php';
		// getting the current size of the user folder in Gigabytes
		$currentSize = getUserUsageSize($user);
		// getting user quota in Gigabytes
		$quota = getUserQuota($user);
		// Setting boolean variable that will indicate whether the user has exceeded it's allocated space, if true the button to add new dataset will not appear
		$exceededSpace = $quota > $currentSize ? FALSE : TRUE;
		if ($exceededSpace)
			echo "<span style='color:#FF0000; font-weight: bold;'>You have exceeded your quota (" . $quota . "G) please clear space and then reload to add new dataset</span><br><br>";

	}
?>
<table width="100%" cellpadding="0"><tr>
<td width="25%" valign="top">
	<?php
	// .------------------.
	// | Make new project |
	// '------------------'
	if (isset($_SESSION['logged_on'])) {
		// show Install new dataset button only if user has space
		if(!$exceededSpace)
			echo "<input name='button_InstallNewDataset' type='button' value='Install New Dataset' onclick='parent.show_hidden(\"Hidden_InstallNewDataset\")'><br>";

		$_SESSION['pending_install_project_count'] = 0;
		?>
		<b><font size='2'>Datasets Pending</font></b><br>
		<div class='tab' style='color:#CC0000; font-size:10pt;' id='newly_installed_list' name='newly_installed_list'></div><br>
		<div style='color:#CC0000; font-size:10pt; visibility:hidden; text-align:left;' id='pending_comment' name='pending_comment'>Reload page after any current uploads have completed to prepare pending datasets
		for upload.<br><br>Additional datasets can be defined with the 'Install New Dataset' button while files upload.<br></div>
		<div style='color:#CC0000; font-size:10pt; visibility:hidden; text-align:center;' id='name_error_comment' name='name_error_comment'>(Entered dataset name is already in use.)</div>
		<?php
	}
	?>
</td>
<td width="75%" valign="top">
	<?php
	// .---------------.
	// | User projects |
	// '---------------'
	$userProjectCount = 0;
	if (isset($_SESSION['logged_on'])) {
		$projectsDir      = "users/".$user."/projects/";
		$projectFolders   = array_diff(glob($projectsDir."*"), array('..', '.'));
		// Sort directories by date, newest first.
		array_multisort(array_map('filemtime', $projectFolders), SORT_DESC, $projectFolders);
		// Trim path from each folder string.
		foreach($projectFolders as $key=>$folder) {   $projectFolders[$key] = str_replace($projectsDir,"",$folder);   }
		// Split project list into ready/working/starting lists for sequential display.
		$projectFolders_complete = array();
		$projectFolders_working  = array();
		$projectFolders_starting = array();
		foreach($projectFolders as $key=>$project) {
			if (file_exists("users/".$user."/projects/".$project."/complete.txt")) {
				array_push($projectFolders_complete,$project);
			} else if (file_exists("users/".$user."/projects/".$project."/working.txt")) {
				array_push($projectFolders_working, $project);
			} else if (is_dir("users/".$user."/projects/".$project)) {
				array_push($projectFolders_starting,$project);
			}
		}
		$userProjectCount_starting = count($projectFolders_starting);
		$userProjectCount_working  = count($projectFolders_working);
		$userProjectCount_complete = count($projectFolders_complete);
		// Sort complete and working projects alphabetically.
		array_multisort($projectFolders_working,  SORT_ASC, $projectFolders_working);
		array_multisort($projectFolders_complete, SORT_ASC, $projectFolders_complete);
		// Build new 'projectFolders' array;
		$projectFolders   = array();
		$projectFolders   = array_merge($projectFolders_starting, $projectFolders_working, $projectFolders_complete);
		$userProjectCount = count($projectFolders);
		// displaying size if it's bigger then 0
		if ($currentSize > 0)
			echo "<b><font size='2'>User installed datasets: (currently using " . $currentSize . "G of " . $quota . "G)</font></b>\n\t\t\t\t";
		else
			echo "<b><font size='2'>User installed datasets:</font></b>\n\t\t\t\t";
		echo "<br>\n\t\t\t\t";
		foreach($projectFolders_starting as $key_=>$project) {
			// Load colors for project.
			$colorFile        = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colorFile)) {
				$handle       = fopen($colorFile,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}
			$parentFile   = "users/".$user."/projects/".$project."/parent.txt";
			$handle       = fopen($parentFile,'r');
			$parentString = trim(fgets($handle));
			fclose($handle);
			// getting project name
			$projectNameString = file_get_contents("users/".$user."/projects/".$project."/name.txt");
			$projectNameString = trim($projectNameString);

			$key = $key_;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#CC0000;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";

			echo $projectNameString;
			echo "</font></span>\n\t\t\t\t";
			echo "<span id='p_delete_".$key."'></span><br>\n\t\t\t\t";
			echo "<div id='frameContainer.p3_".$key."'></div>\n";
		}
		foreach($projectFolders_working as $key_=>$project) {
			// Load colors for project.
			$colorFile        = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colorFile)) {
				$handle       = fopen($colorFile,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}
			$parentFile   = "users/".$user."/projects/".$project."/parent.txt";
			$handle       = fopen($parentFile,'r');
			$parentString = trim(fgets($handle));
			fclose($handle);
			// getting project name
			$projectNameString = file_get_contents("users/".$user."/projects/".$project."/name.txt");
			$projectNameString = trim($projectNameString);

			$key = $key_ + $userProjectCount_starting;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#BB9900;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";
			echo $projectNameString;

			echo "</font></span>\n\t\t\t\t";
			echo "<span id='p_delete_".$key."'></span><br>\n\t\t\t\t";
			echo "<div id='frameContainer.p2_".$key."'></div>\n";
		}
?>
<script type='text/javascript'>
	function loadExternal(imageUrl) {
		window.open(imageUrl);
	}
</script>
<?php
		foreach($projectFolders_complete as $key_=>$project) {
			// Load data type for project.
			$handle     = fopen("users/".$user."/projects/".$project."/dataFormat.txt", "r");
			$dataFormat = fgets($handle);
			fclose($handle);

			// Load colors for project.
			$colorFile        = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colorFile)) {
				$handle       = fopen($colorFile,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}
			$parentFile   = "users/".$user."/projects/".$project."/parent.txt";
			$handle       = fopen($parentFile,'r');
			$parentString = trim(fgets($handle));
			fclose($handle);
			// getting project name
			$projectNameString = file_get_contents("users/".$user."/projects/".$project."/name.txt");
			$projectNameString = trim($projectNameString);

			$key = $key_ + $userProjectCount_starting + $userProjectCount_working;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#00AA00;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).". ";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";
			echo $projectNameString."</font>";

			// display total project size
			$totalSizeFile = "users/".$user."/projects/".$project. "/totalSize.txt";
			// first checking if size already calculated and is stored in totalSize.txt
			if (file_exists($totalSizeFile))
			{
				$handle       = fopen($totalSizeFile,'r');
				$projectSizeStr = trim(fgets($handle));
				fclose($handle);
			}
			else // calculate size and store in totalSize.txt to avoid calculating again
			{
				// calculating size
				$projectSizeStr = trim(shell_exec("du -sh " . "users/".$user."/projects/".$project. "/ | cut -f1"));
				// saving to file
				$output       = fopen($totalSizeFile, 'w');
				fwrite($output, $projectSizeStr);
				fclose($output);
			}
			// printing total size
			echo " <font color='black' size='1'>(". $projectSizeStr .")</font>";

			echo "</span>";
			echo "<span id='p_delete_".$key."'></span><br>\n\t\t\t\t";
			echo "<div id='frameContainer.p1_".$key."'></div>\n";
		}
	}
	?>
</td>
</tr></table>

<?php
	//.-----------------.
	//| System projects |
	//'-----------------'
	$projectsDir          = "users/default/projects/";
	$systemProjectFolders = array_diff(glob($projectsDir."*"), array('..', '.'));
	// Sort directories by date, newest first.
	array_multisort(array_map('filemtime', $systemProjectFolders), SORT_DESC, $systemProjectFolders);
	// Trim path from each folder string.
	foreach($systemProjectFolders as $key=>$folder) {   $systemProjectFolders[$key] = str_replace($projectsDir,"",$folder);   }
	$systemProjectCount = count($systemProjectFolders);
?>


<script type="text/javascript">
var userProjectCount   = "<?php echo $userProjectCount; ?>";
var systemProjectCount = "<?php echo $systemProjectCount; ?>";
<?php
//.----------------.
//| Uploading data |
//'----------------'
if (isset($_SESSION['logged_on'])) {
	foreach($projectFolders_starting as $key_=>$project) {    // frameContainer.p3_[$key] : starting.
		$key      = $key_;
		$project  = $projectFolders[$key];
		// Read in dataFormat string for project.
		$handle     = fopen("users/".$user."/projects/".$project."/dataFormat.txt", "r");
		$dataFormat = fgets($handle);
		fclose($handle);
		echo "\n// javascript for project #".$key.", '".$project."'\n";
		echo "var el_p               = document.getElementById('frameContainer.p3_".$key."');\n";
		// Javascript to build file load button interface.
		echo "el_p.innerHTML         = '<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" ";
		if ((strlen($dataFormat) > 1) && ($dataFormat[2] == '1')) {
			// paired files to be uploaded.
			echo "style=\"height:76px\" src=\"uploader.2.php\"";
		} else {
			// single file to be uploaded.
			echo "style=\"height:38px\" src=\"uploader.1.php\"";
		}
		echo " marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\" width=\"100%\" frameborder=\"0\"></iframe>';\n";
		echo "var p_iframe           = document.getElementById('p_".$key."');\n";
		echo "var p_js               = p_iframe.contentWindow;\n";
		echo "p_js.display_string    = new Array();\n";
		echo "p_js.user              = '".$user."';\n";
		echo "p_js.project           = '".$project."';\n";
		echo "p_js.key               = 'p_".$key."';\n";
		if ($dataFormat == '0') {
			// SnpCgh microarray
			echo "p_js.display_string[0] = 'Add : SnpCgh array data...';\n";
			echo "p_js.dataFormat        = 'SnpCghArray';\n";
		} else if (($dataFormat == '1:0:0') || ($dataFormat == '1:0:1')) {
			// WGseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read WGseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'WGseq_single';\n";
		} else if (($dataFormat == '1:1:0') || ($dataFormat == '1:1:1')) {
			// WGseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read WGseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read WGseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'WGseq_paired';\n";
		} else if (($dataFormat == '1:2:0') || ($dataFormat == '1:2:1') || ($dataFormat == '1:3:0') || ($dataFormat == '1:3:1')) {
			// WGseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : WGseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.dataFormat        = 'WGseq_single';\n";
		} else if (($dataFormat == '2:0:0') || ($dataFormat == '2:0:1')) {
			// ddRADseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read ddRADseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'ddRADseq_single';\n";
		} else if (($dataFormat == '2:1:0') || ($dataFormat == '2:1:1')) {
			// ddRADseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read ddRADseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read ddRADseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'ddRADseq_paired';\n";
		} else if (($dataFormat == '2:2:0') || ($dataFormat == '2:2:1') || ($dataFormat == '2:3:0') || ($dataFormat == '2:3:1')) {
			// ddRADseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : ddRADseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.dataFormat        = 'ddRADseq_single';\n";
		} else if (($dataFormat == '3:0:0') || ($dataFormat == '3:0:1')) {
			// RNAseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read RNAseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'RNAseq_single';\n";
		} else if (($dataFormat == '3:1:0') || ($dataFormat == '3:1:1')) {
			// RNAseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read RNAseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read RNAseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'RNAseq_paired';\n";
		} else if (($dataFormat == '3:2:0') || ($dataFormat == '3:2:1') || ($dataFormat == '3:3:0') || ($dataFormat == '3:3:1')) {
			// RNAseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : RNAseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.dataFormat        = 'RNAseq_single';\n";
		} else if (($dataFormat == '4:0:0') || ($dataFormat == '4:0:1')) {
			// IonExpressSeq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read IonExpress data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'IonExpressSeq_single';\n";
		} else if (($dataFormat == '4:1:0') || ($dataFormat == '4:1:1')) {
			// IonExpressSeq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read IonExpress data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read IonExpress data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.dataFormat        = 'IonExpressSeq_paired';\n";
		} else if (($dataFormat == '4:2:0') || ($dataFormat == '4:2:1') || ($dataFormat == '4:3:0') || ($dataFormat == '4:3:1')) {
			// IonExpressSeq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : IonExpress data (SAM/BAM/TXT)...';\n";
			echo "p_js.dataFormat        = 'IonExpressSeq_single';\n";
		}
	}
	foreach($projectFolders_working as $key_=>$project) {   // frameContainer.p2_[$key] : working.
		$key      = $key_ + $userProjectCount_starting;
		$project  = $projectFolders[$key];
		$handle   = fopen("users/".$user."/projects/".$project."/dataFormat.txt", "r");
		$dataFormat = fgets($handle);
		fclose($handle);
		echo "\n// javascript for project #".$key.", '".$project."'\n";
		echo "var el_p            = document.getElementById('frameContainer.p2_".$key."');\n";
		echo "el_p.innerHTML      = '<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px; border:0px;\" ";
		echo     "src=\"project.working.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\" width=\"100%\" frameborder=\"0\"></iframe>';\n";
		echo "var p_iframe        = document.getElementById('p_".$key."');\n\t";
		echo "var p_js            = p_iframe.contentWindow;\n";
		echo "p_js.user           = \"".$user."\";\n\t";
		echo "p_js.project        = \"".$project."\";\n\t";
		echo "p_js.key            = \"p_".$key."\";\n";
	}
}
?>
</script>
