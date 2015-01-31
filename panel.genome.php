<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>Install new reference genomes for use in sequence analysis.</font><br><br>
<table width="100%" cellpadding="0"><tr>
<td width="50%">
	<?php
	// .-----------------.
	// | Make new genome |
	// '-----------------'
	if (isset($_SESSION['logged_on'])) {
		echo "<input name=\"button_InstallNewGenome\"  type=\"button\" value=\"Install New Genome\"  onclick=\"parent.show_hidden('Hidden_InstallNewGenome')\"><br>\n\t\t\t\t";

		$_SESSION['pending_install_genome_count'] = 0;
		?>
		<b><font size='2'>Genomes Pending</font></b><br>
		<div class='tab' style='color:#CC0000; font-size:10pt;' id='newly_installed_list' name='newly_installed_list'></div><br>
		<div style='color:#CC0000; font-size:10pt; visibility:hidden; text-align:center;' id='pending_comment' name='pending_comment'>(Reload page after current uploads have completed to prepare these
		for upload.)</div>
		<?php
	}

	// .--------------.
	// | User genomes |
	// '--------------'
	$userGenomeCount = 0;
	if (isset($_SESSION['logged_on'])) {
		$genomesDir    = "users/".$user."/genomes/";
		$genomeFolders = array_diff(glob($genomesDir."*"), array('..', '.'));
		// Sort directories by date, newest first.
		array_multisort(array_map('filemtime', $genomeFolders), SORT_DESC, $genomeFolders);
		// Trim path from each folder string.
		foreach($genomeFolders as $key=>$folder) {   $genomeFolders[$key] = str_replace($genomesDir,"",$folder);   }
		// Split genome list into ready/working/starting lists for sequential display.
		$genomeFolders_complete = array();
		$genomeFolders_working  = array();
		$genomeFolders_starting = array();
		foreach($genomeFolders as $key=>$genome) {
			if (file_exists("users/".$user."/genomes/".$genome."/complete.txt")) {
				array_push($genomeFolders_complete,$genome);
			} else if (file_exists("users/".$user."/genomes/".$genome."/working.txt")) {
				array_push($genomeFolders_working, $genome);
			} else {
				array_push($genomeFolders_starting,$genome);
			}
		}
		$userGenomeCount_starting = count($genomeFolders_starting);
		$userGenomeCount_working  = count($genomeFolders_working);
		$userGenomeCount_complete = count($genomeFolders_complete);
		// Sort complete and working genomes alphabetically.
		array_multisort($genomeFolders_working,  SORT_ASC, $genomeFolders_working);
		array_multisort($genomeFolders_complete, SORT_ASC, $genomeFolders_complete);
		// Build new 'genomeFolders' array;
		$genomeFolders   = array();
		$genomeFolders   = array_merge($genomeFolders_starting, $genomeFolders_working, $genomeFolders_complete);
		$userGenomeCount = count($genomeFolders);

		echo "<b><font size='2'>User Installed Genomes:</font></b>\n\t\t\t\t";
		echo "<br>\n\t\t\t\t";
		foreach($genomeFolders_starting as $key_=>$genome) {
			$genomeNameString = file_get_contents("users/".$user."/genomes/".$genome."/name.txt");
			$genomeNameString = trim($genomeNameString);
			$key = $key_;
			echo "<span id='g_label_".$key."' style='color:#CC0000;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='genome_delete_".$key."' type='button' onclick=\"parent.deleteGenomeConfirmation('".$user."','".$genome."','".$key."')\">Delete</button>";
			echo $genomeNameString;

			$sizeFile_1   = "users/".$user."/genomes/".$genome."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='g_size1_".$key."'></span>"; }

			echo "</font></span>\n\t\t\t\t";
			echo "<span id='g_delete_".$key."'></span>\n\t\t";
			echo "\n\t\t\t\t";
			echo "<div id='frameContainer.g3_".$key."'></div>";
		}
		foreach($genomeFolders_working as $key_=>$genome) {
			$genomeNameString = file_get_contents("users/".$user."/genomes/".$genome."/name.txt");
			$genomeNameString = trim($genomeNameString);
			$key = $key_ + $userGenomeCount_starting;
			echo "<span id='g_label_".$key."' style='color:#BB9900;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<button id='genome_delete_".$key."' type='button' onclick=\"parent.deleteGenomeConfirmation('".$user."','".$genome."','".$key."')\">Delete</button>";
			echo $genomeNameString;

			$sizeFile_1   = "users/".$user."/genomes/".$genome."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='g_size1_".$key."'></span>"; }

			echo "</font></span>\n\t\t\t\t";
			echo "<span id='g_delete_".$key."'></span>\n\t\t";
			echo "\n\t\t\t\t";
			echo "<div id='frameContainer.g2_".$key."'></div>";
		}
		foreach($genomeFolders_complete as $key_=>$genome) {
			$genomeNameString = file_get_contents("users/".$user."/genomes/".$genome."/name.txt");
			$genomeNameString = trim($genomeNameString);
			$key = $key_ + $userGenomeCount_starting + $userGenomeCount_working;
			echo "<span id='g_label_".$key."' style='color:#00AA00;'>\n\t\t\t\t";
			echo "<font size='2'>".($key+1).". ";
			echo "<button id='genome_delete_".$key."' type='button' onclick=\"parent.deleteGenomeConfirmation('".$user."','".$genome."','".$key."')\">Delete</button>";
			echo $genomeNameString;

			$sizeFile_1   = "users/".$user."/genomes/".$genome."/upload_size_1.txt";
			$handle       = fopen($sizeFile_1,'r');
			$sizeString_1 = trim(fgets($handle));
			fclose($handle);
			if ($sizeString_1 !== "") { echo " <font color='black' size='1'>(".$sizeString_1." bytes)</font>";
			} else {                    echo " <span id='g_size1_".$key."'></span>"; }

			echo "</font></span>\n\t\t\t\t";
			echo "<span id='g_delete_".$key."'></span>\n\t\t";
			echo "\n\t\t\t\t";
			echo "<div id='frameContainer.g1_".$key."'></div>";
		}
	}
	?>
</td><td width="50%" valign="top">
	<?php
	//.----------------.
	//| System genomes |
	//'----------------'
	$genomesDir          = "users/default/genomes/";
	$systemGenomeFolders = array_diff(glob($genomesDir."*"), array('..', '.'));
	// Sort directories by date, newest first.
	array_multisort(array_map('filemtime', $systemGenomeFolders), SORT_DESC, $systemGenomeFolders);
	// Trim path from each folder string.
	foreach($systemGenomeFolders as $key=>$folder) {   $systemGenomeFolders[$key] = str_replace($genomesDir,"",$folder);   }
	$systemGenomeCount = count($systemGenomeFolders);
	echo "<b><font size='2'>Installed Reference Genomes:</font></b>\n";
	foreach ($systemGenomeFolders as $key=>$genome) {
		$genomeNameString = file_get_contents("users/default/genomes/".$genome."/name.txt");
		$genomeNameString = trim($genomeNameString);
		echo "\t\t\t\t<font size='2'>".$genomeNameString."</font><br>\n";
	}
	?>
</td></tr></table>

<script type="text/javascript">
// Javascript which needs to be run.
var userGenomeCount   = "<?php echo $userGenomeCount; ?>";
var systemGenomeCount = "<?php echo $systemGenomeCount; ?>";
<?php
if (isset($_SESSION['logged_on'])) {
	foreach($genomeFolders_starting as $key_=>$genome) {
		$key      = $key_;
		$genome   = $genomeFolders[$key];
		echo "\n\t// javascript for genome #".$key.", '".$genome."'\n\t";
		echo "var el_g               = document.getElementById('frameContainer.g3_".$key."');\n\t";
		echo "el_g.innerHTML         = '<iframe id=\"g_".$key."\" name=\"g_".$key."\" class=\"upload\" style=\"height:28px\" ";
		echo "src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\" width=\"100%\" frameborder=\"0\"></iframe>';\n\t";
		echo "var g_iframe           = document.getElementById('g_".$key."');\n\t";
		echo "var g_js               = g_iframe.contentWindow;\n\t";
		echo "g_js.display_string    = new Array();\n\t";
		echo "g_js.display_string[0] = \"Add : Genome reference FASTA file...\";\n\t";
		echo "g_js.target_dir        = \"users/".$user."/genomes/".$genome."/\";\n\t";
		echo "g_js.conclusion_script = \"scripts_genomes/genome.install_1.php\";\n\t";
		echo "g_js.user              = \"".$user."\";\n\t";
		echo "g_js.genome            = \"".$genome."\";\n\t";
		echo "g_js.key               = \"g_".$key."\";\n";
	}
	foreach($genomeFolders_working as $key_=>$genome) {
		$key      = $key_ + $userGenomeCount_starting;
		$genome   = $genomeFolders[$key];
		echo "\n\t// javascript for genome #".$key.", '".$genome."'\n\t";
		echo "var el_g               = document.getElementById('frameContainer.g2_".$key."');\n\t";
		echo "el_g.innerHTML         = '<iframe id=\"g_".$key."\" name=\"g_".$key."\" class=\"upload\" style=\"height:76px\" ";
		echo "src=\"genome.working.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\" width=\"100%\" frameborder=\"0\"></iframe>';\n\t";
		echo "var g_iframe           = document.getElementById('g_".$key."');\n\t";
		echo "var g_js               = g_iframe.contentWindow;\n\t";
		echo "g_js.user              = \"".$user."\";\n\t";
		echo "g_js.genome            = \"".$genome."\";\n\t";
		echo "g_js.key               = \"g_".$key."\";\n";
	}
}
?>
</script>
