$token = "$Env:av_token"
$headers = @{}
$headers['Authorization'] = "Bearer $token"
$headers["Content-type"] = "application/json"
$body = @{}
$body['accountName'] = "rickhelmus"
$body['projectSlug'] = "patRoonDeps"
$bodyAsJson = $body | ConvertTo-json
Invoke-RestMethod -Uri 'https://ci.appveyor.com/api/account/rickhelmus/builds' -Headers $headers -Method Post -Body $bodyAsJson -ContentType "application/json"