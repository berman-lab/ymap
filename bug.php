<?php
    session_start();
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
        "http://www.w3.org/TR/html4/loose.dtd">
<?php
	require_once 'php/constants.php';
	chmod($GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'], 0777);
?>
<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Yeast Haplotype Analysis System</title>
		<script type="text/javascript" src="js/jquery-1.7.2.js"></script>
		<script type="text/javascript" src="js/jquery.form.js"></script>
		<script>
				function submitBug(project){
					submitter = $('#submitter').val();
					description = $('#description').val();
					$.ajax({
						url : '<?php echo $GLOBALS['url']?>php/addBug_server.php',
						type : 'post',
						data : {
							submitter: submitter,
							description: description
						},
						success : function(answer){
							window.location.href=window.location.href;
						}
					});
				}
		</script>
	</head>
	<body>
		<div id="bugSubmit">
			<button type="button" onclick="window.location.href='index.php'">Back to Home</button><b>User: <?php echo $_SESSION['user']?></b>
			<input type="hidden" id="submitter" name="submitter" value="<?php echo $_SESSION['user']?>"><br>
			<textarea id="description" rows="15" cols="150">Describe your bug or feature request.</textarea><br />
			<button type="button" onclick="submitBug()">Submit Bug or Feature!</button><br /><br />
		</div>
		<div id="bugList">
			<table id="bugTable" border="1">
				<tr>
					<th>Submitter</th>
					<th colspan=80>Description</th>
					<th>Status</th>
					<th colspan=80>Notes</th>
				</tr>
					<?php
					$bugtrackFile = $GLOBALS['directory'].$GLOBALS['bugtrackDir'].$GLOBALS['bugtrackFile'];
					if (file_exists($bugtrackFile)){
						$bugFileHandle = fopen($bugtrackFile, 'r');
						while (($buffer = fgets($bugFileHandle, 4096)) !== false) {
							list ($submitter, $description, $status, $notes) = explode("|", $buffer);
							echo "
								<tr>
									<td>".$submitter."</td>
									<td colspan=80>".$description."</td>
									<td>".$status."</td>
									<td colspan=80>".$notes."</td>
								</tr>
							";
						}
						fclose($bugFileHandle);
					} else {
						$bugsOutput = fopen($bugtrackFile, 'w');
						chmod($bugsOutput, 0777);
						fwrite($bugtrackFile, "test|test|test|test|test");
						fclose($bugsOutput);
					}
					?>
			</table>
		</div>
	</body>
</html>
