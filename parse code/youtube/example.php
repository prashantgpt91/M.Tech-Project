<?php
  $q = "tseries";
  $maxResults = 50;
  require_once ('google-api-php-client/src/Google_Client.php');
  require_once ('google-api-php-client/src/contrib/Google_YouTubeService.php');

  $DEVELOPER_KEY = 'AIzaSyC-qZc4ERxavJ_d0rpnxwSekjklc0Z1slY';

  $client = new Google_Client();
  $client->setDeveloperKey($DEVELOPER_KEY);

  $youtube = new Google_YoutubeService($client);

 
  
    $searchResponse = $youtube->search->listSearch('id,snippet', array(
      'q' => $q,
      'maxResults' => $maxResults,
    ));
	 echo json_decode($searchRespons);
    $videos = '';
    $channels = '';

    foreach ($searchResponse['items'] as $searchResult) {
	echo json_encode($searchResult)."<br><br><br><br><br><br><br>";
      	echo $searchResult['snippet']['channelTitle']."<br>";
      switch ($searchResult['id']['kind']) {
        case 'youtube#video':
          $videos .= sprintf('<li>%s (%s)</li>', $searchResult['snippet']['channelTitle'],
            $searchResult['id']['videoId']."<a href=http://www.youtube.com/watch?v=".$searchResult['id']['videoId']." target=_blank>   Watch This Video</a>");
          break;
        case 'youtube#channel':
          $channels .= sprintf('<li>%s (%s)</li>', $searchResult['snippet']['title'],
            $searchResult['id']['channelId']);
          break;
       }
    }

   

?>