<?php
    session_start();
    require_once 'constants.php';

//print_r($_SERVER);
//print_r($GLOBALS);

$file_path = $_SERVER[SCRIPT_FILENAME];
$file_dir  = dirname($file_path);
$file_name = basename($file_path);

$files = glob($file_dir.'/files/*'); // get all file names
foreach($files as $file){ // iterate files
  if(is_file($file))
    unlink($file); // delete file
}

$files = glob($file_dir.'/files2/*'); // get all file names
foreach($files as $file){ // iterate files
  if(is_file($file))
    unlink($file); // delete file
}

header("Location: ../resetter.html");
exit;

?>
