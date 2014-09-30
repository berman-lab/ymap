<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>Install SNP/CGH microarray and next generation sequence datasets to, or delete datasets from, your user account.</font><br><br>
<table width="100%" cellpadding="0"><tr>
<td width="25%" valign="top">
	<?php
	// .------------------.
	// | Make new project |
	// '------------------'
	if (isset($_SESSION['logged_on'])) {
		echo "<input name='button_InstallNewDataset' type='button' value='Install New Dataset' onclick='parent.show_hidden(\"Hidden_InstallNewDataset\")'><br>";

		$_SESSION['pending_install_project_count'] = 0;
		?>
		<b><font size='2'>Datasets Pending</font></b><br>
		<div class='tab' style='color:#CC0000; font-size:10pt;' id='newly_installed_list' name='newly_installed_list'></div><br>
		<div style='color:#CC0000; font-size:10pt; visibility:hidden; text-align:left;' id='pending_comment' name='pending_comment'>Reload page after any current uploads have completed to prepare pending datasets
		for upload.<br><br>Additional datasets can be defined with the 'Install New Dataset' button while files upload.</div>
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
			} else {
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

			$key = $key_;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#CC0000;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";
			echo $project;

			$sizeFile_1   = "users/".$user."/projects/".$project."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			$sizeFile_2   = "users/".$user."/projects/".$project."/upload_size_2.txt";
			$handle       = fopen($sizeFile_2,'r');
			$sizeString_2 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='p_size1_".$key."'></span>"; }
			if ($sizeString_2 !== "") { echo " <font color='black' size='1'>(".$sizeString_2." bytes)</font>";
			} else {                    echo " <span id='p_size2_".$key."'></span>"; }

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

			$key = $key_ + $userProjectCount_starting;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#BB9900;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";
			echo $project;

			$sizeFile_1   = "users/".$user."/projects/".$project."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			$sizeFile_2   = "users/".$user."/projects/".$project."/upload_size_2.txt";
			$handle       = fopen($sizeFile_2,'r');
			$sizeString_2 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='p_size1_".$key."'></span>"; }
			if ($sizeString_2 !== "") { echo " <font color='black' size='1'>(".$sizeString_2." bytes)</font>";
			} else {                    echo " <span id='p_size2_".$key."'></span>"; }

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

			$key = $key_ + $userProjectCount_starting + $userProjectCount_working;
			echo "\n\t\t\t\t<!-- project '{$project}', #{$key}. --!>\n\t\t\t\t";
			echo "<span id='p_label_".$key."' style='color:#00AA00;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).". ";
			echo "<button id='project_delete_".$key."' type='button' onclick=\"parent.deleteProjectConfirmation('".$user."','".$project."','".$key."')\">Delete</button>";
			echo $project."</font>";

			$sizeFile_1   = "users/".$user."/projects/".$project."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			$sizeFile_2   = "users/".$user."/projects/".$project."/upload_size_2.txt";
			$handle       = fopen($sizeFile_2,'r');
			$sizeString_2 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='p_size1_".$key."'></span>"; }
			if ($sizeString_2 !== "") { echo " <font color='black' size='1'>(".$sizeString_2." bytes)</font>";
			} else {                    echo " <span id='p_size2_".$key."'></span>"; }

			echo "</span> : ";
			echo "<span onclick='loadExternal(\"users/".$user."/projects/".$project."/SNP_CNV_v1.txt\")'><font size='1'>[SNP/CNV data]</font></span> ";
			echo "<span onclick='loadExternal(\"users/".$user."/projects/".$project."/putative_SNPs_v4.txt\")'><font size='1'>[SNP data]</font></span>\n\t\t\t\t";
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
if (isset($_SESSION['logged_on'])) {
	foreach($projectFolders_starting as $key_=>$project) {    // frameContainer.p3_[$key] : starting.
		$key      = $key_;
		$project  = $projectFolders[$key];
		$handle   = fopen("users/".$user."/projects/".$project."/dataType.txt", "r");
		$dataType = fgets($handle);
		fclose($handle);
		echo "\n// javascript for project #".$key.", '".$project."'\n";
		echo "var el_p            = document.getElementById('frameContainer.p3_".$key."');\n";
//		echo "console.log(el_p);";
		echo "el_p.innerHTML      = '<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" ";
		if ((strlen($dataType) > 1) && ($dataType[2] == '1')) {   echo "style=\"height:76px\" src=\"uploader.2.php\"";   } else {   echo "style=\"height:38px\" src=\"uploader.1.php\"";   }
		echo " marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\" width=\"100%\" frameborder=\"0\"></iframe>';\n";
		echo "var p_iframe        = document.getElementById('p_".$key."');\n";
		echo "var p_js            = p_iframe.contentWindow;\n";
		echo "p_js.display_string = new Array();\n";
		echo "p_js.target_dir     = '../../users/".$user."/projects/".$project."/';\n";
		echo "p_js.user           = '".$user."';\n";
		echo "p_js.project        = '".$project."';\n";
		echo "p_js.key            = 'p_".$key."';\n";
		if ($dataType == '0') {                                      // SnpCgh microarray
			echo "p_js.display_string[0] = 'Add : SnpCgh array data...';\n";
			echo "p_js.conclusion_script = 'php/project.SnpCgh.install.php';\n";
		} else if ($dataType == '1:0') {                             // WGseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read WGseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_WGseq.install_1.php';\n";
		} else if ($dataType == '1:1') {                             // WGseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read WGseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read WGseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.paired_WGseq.install_1.php';\n";
		} else if (($dataType == '1:2') || ($dataType == '1:3')) {   // WGseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : WGseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_WGseq.install_1.php';\n";
		} else if ($dataType == '2:0') {                             // ddRADseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read ddRADseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_ddRADseq.install_1.php';\n";
		} else if ($dataType == '2:1') {                             // ddRADseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read ddRADseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read ddRADseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.paired_ddRADseq.install_1.php';\n";
		} else if (($dataType == '2:2') || ($dataType == '2:3')) {   // ddRADseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : ddRADseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_ddRADseq.install_1.php';\n";
		} else if ($dataType == '3:0') {                             // RNAseq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read RNAseq data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_RNAseq.install_1.php';\n";
		} else if ($dataType == '3:1') {                             // RNAseq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read RNAseq data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read RNAseq data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.paired_RNAseq.install_1.php';\n";
		} else if (($dataType == '3:2') || ($dataType == '3:3')) {   // RNAseq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : RNAseq data (SAM/BAM/TXT)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_RNAseq.install_1.php';\n";
		} else if ($dataType == '4:0') {                             // IonExpressSeq : single-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Single-end-read IonExpress data (FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_IonExpressSeq.install_1.php';\n";
		} else if ($dataType == '4:1') {                             // IonExpressSeq : paired-end [FASTQ/ZIP/GZ]
			echo "p_js.display_string[0] = 'Add : Paired-end-read IonExpress data (1/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.display_string[1] = 'Add : Paired-end-read IonExpress data (2/2; FASTQ/ZIP/GZ)...';\n";
			echo "p_js.conclusion_script = 'php/project.paired_IonExpressSeq.install_1.php';\n";
		} else if (($dataType == '4:2') || ($dataType == '4:3')) {   // IonExpressSeq : [SAM/BAM/TXT]
			echo "p_js.display_string[0] = 'Add : IonExpress data (SAM/BAM/TXT)...';\n";
			echo "p_js.conclusion_script = 'php/project.single_IonExpressSeq.install_1.php';\n";
		}
	}
	foreach($projectFolders_working as $key_=>$project) {   // frameContainer.p2_[$key] : working.
		$key      = $key_ + $userProjectCount_starting;
		$project  = $projectFolders[$key];
		$handle   = fopen("users/".$user."/projects/".$project."/dataType.txt", "r");
		$dataType = fgets($handle);
		fclose($handle);
		echo "\n// javascript for project #".$key.", '".$project."'\n";
		echo "var el_p            = document.getElementById('frameContainer.p2_".$key."');\n";
//		echo "console.log(el_p);";
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
