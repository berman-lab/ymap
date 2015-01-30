<?php
    session_start();
    if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
    require_once 'constants.php';
    echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
?>
<html lang="en">
	<head>
		<style type="text/css">
			body {font-family: arial;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>[Needs Title]</title>
	</head>
	<body>
		<div id="genomeCreationInformation"><p>
			<form action="genome.create_server.test.php" method="post">
				<label for="newGenomeName">Genome Name: </label><input type="text" name="newGenomeName" id="newGenomeName" size="100"><br>
				<br>
				<input type="submit" value="Create New Genome">
			</form>
		</p></div>
	</body>
</html>
