<?php 
session_start();
require_once 'constants.php';

session_destroy();
header('Location: '.$GLOBALS['url']);
?>
