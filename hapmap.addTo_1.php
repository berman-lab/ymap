<?php
	session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">
<?php
	//session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: user.login.php');
	}

	$bad_chars = array(".", ",", "\\", "/", " ");
	$hapmap    = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$user      = filter_input(INPUT_POST, "user",   FILTER_SANITIZE_STRING);
	$key       = filter_input(INPUT_POST, "key",    FILTER_SANITIZE_STRING);

	$dir1      = "users/".$user."/hapmaps";
	$dir2      = "users/".$user."/hapmaps/".$hapmap;
	$dir3      = "users/default/hapmaps/".$hapmap;

	// figure out what user the hapmap is installed under.
	$folder = "users/".$user."/hapmaps/".$hapmap."/";

	// Load genome from 'hapmap/genome.txt'.
	$handle1 = fopen($folder."genome.txt", "r");
	$genome  = trim(fgets($handle1));
	fclose($handle1);

	// Load parent from 'hapmap/parent.txt'.
	$handle2 = fopen($folder."parent.txt", "r");
	$parent  = trim(fgets($handle2));
	fclose($handle2);
?>
<html lang="en">
	<HEAD>
		<style type="text/css">
			body {font-family: arial;}
			.tab {margin-left:   1cm;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>[Needs Title]</title>
	</HEAD>
	<BODY onload="UpdateNextList()">
		<div id="loginControls"><p>
			<?php
			if (isset($_SESSION['logged_on'])) {
				$user = $_SESSION['user'];
				echo "user: ".$user."<br>";
			}
			?>

			<!-- If user is logged in, show logout button, otherwise, show the login button so we can get the user logged in-->
			<button type="button" onclick="window.location.href='<?php if(isset($_SESSION['logged_on'])){echo "logout_server.php";}else{echo "user.login.php";}?>'"><?php if(isset($_SESSION['logged_on'])){echo "Logout";}else{echo "Login";}?></button>
			<button type="button" onclick="window.location.href='index.php'">Back to Home</button>
		</p></div>
		<div id="hapmapCreationInformation"><p>
			<form action="hapmap.addTo_2.php" method="post">
				<table><tr bgcolor="#CCFFCC"><td>
					<label for="hapmap">Hapmap Name : </label><input type="text" name="hapmap" id="hapmap" value="<?php echo $hapmap; ?>" readonly>
				</td><td>
					Unique name for this hapmap.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="genome">Reference genome : </label><input type="text" name="genome" id="genome" value="<?php echo $genome; ?>" size="60" readonly>
				</td><td valign="top">
					Reference genome used to construct hapmap.
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<label for="parent">Parental strain : </label><input type="text" name="parent" id="parent" value="<?php echo $parent; ?>" readonly>
				</td><td valign="top">
					Only whole genome sequence datasets are used in constructing hapmaps.
				</td></tr><tr bgcolor="#CCCCFF"><td valign="top">
					<?php
					// figure out which hapmaps have been defined for this species, if any.
					$projectsDir1       = "users/default/projects/";
					$projectsDir2       = "users/".$user."/projects/";
					$projectFolders1    = array_diff(glob($projectsDir1."*"), array('..', '.'));
					$projectFolders2    = array_diff(glob($projectsDir2."*"), array('..', '.'));
					$projectFolders_raw = array_merge($projectFolders1,$projectFolders2);
					// Go through each $projectFolder and look at 'genome.txt' and 'dataType.txt'; build javascript array of prejectName:genome:datatype triplets.
					?>
					Next strain : <select id="selectNext" name="selectNext"><option>[choose]</option></select>
					<script type="text/javascript">
					var nextGenomeDatatype_entries = [['next','genome','dataType']<?php
					foreach ($projectFolders_raw as $key=>$folder) {
						$handle1         = fopen($folder."/genome.txt", "r");
						$genome_string   = trim(fgets($handle1));
						fclose($handle1);
						$handle2         = fopen($folder."/dataType.txt", "r");
						$dataType_string = trim(fgets($handle2));
						$dataType_string = explode(":",$dataType_string);
						$dataType_string = $dataType_string[0];
						fclose($handle2);
						$nextName        = $folder;
						$nextName        = str_replace($projectsDir1,"",$nextName);
						$nextName        = str_replace($projectsDir2,"",$nextName);
						echo ",['{$nextName}','{$genome_string}',{$dataType_string}]";
					}
					?>];

					UpdateNextList=function() {
						var selectedGenome   = "<?php echo $genome; ?>";
						var selectedDatatype = 1
						var select           = document.getElementById("selectNext");     // grab select list.
						select.innerHTML     = '';
						for (var i = 1; i < nextGenomeDatatype_entries.length; i++) {
							var item = nextGenomeDatatype_entries[i];
							if (selectedGenome == item[1] && selectedDatatype == item[2]) {
								var el         = document.createElement("option");
								el.textContent = item[0];
								el.value       = item[0];
								select.appendChild(el);
							}
						}
					}
					</script>
				</td><td valign="top">
					The first dataset used to construct the hapmap. (Others can be added later...)<br>
					Each strain used to construct the hapmap should have large loss of heterozygosity regions.
				</td></tr></table><br>
				<input type="submit" value="Add hapmap entry">
			</form>
		</p></div>
	</body>
</html>
