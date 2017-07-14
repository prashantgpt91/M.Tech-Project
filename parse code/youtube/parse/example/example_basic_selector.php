<?php
// example of how to use basic selector to retrieve HTML contents
include('../simple_html_dom.php');
 
// get DOM from URL or file
$html = file_get_html('1.html');


// find all div tags with id=gbar
//foreach($html->find('li span.yt-channel-title-icon-verified') as $e)
//   echo $e->prev_sibling()->innertext . '<br>';

foreach($html->find('span.yt-channel-title-icon-verified') as $e)
    echo $e->prev_sibling()->innertext."<br>".$e->parent()->parent()->parent()->prev_sibling()->first_child()->innertext . '<br><br>';

?>