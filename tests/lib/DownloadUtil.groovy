
@Grab(group='org.apache.httpcomponents', module='httpclient', version='4.5.13')
@Grab(group='org.apache.httpcomponents', module='httpmime', version='4.5.13')

import org.apache.http.client.methods.HttpGet
import org.apache.http.impl.client.HttpClients
import org.apache.http.client.config.RequestConfig
import org.apache.http.client.utils.URIBuilder
import org.apache.http.entity.ContentType
import org.apache.http.HttpEntity
import org.apache.http.util.EntityUtils

import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.util.zip.ZipInputStream
import java.io.FileOutputStream

public class DownloadUtil {

    public static void downloadAndExtractZip(String url, String destinationFolder) {
        // Create HTTP client
        def httpClient = HttpClients.custom().setDefaultRequestConfig(RequestConfig.custom().setConnectTimeout(5000).build()).build()

        try {
            // Prepare HTTP GET request
            def uri = new URIBuilder(url).build()
            def httpGet = new HttpGet(uri)
            
            // Execute HTTP request
            def response = httpClient.execute(httpGet)
            
            // Check if response is successful
            if (response.getStatusLine().getStatusCode() == 200) {
                // Get HTTP entity
                def entity = response.getEntity()
                
                // Create temporary file to save the downloaded zip
                def tempFile = Files.createTempFile("temp-", ".zip")
                
                // Save the zip file
                Files.copy(entity.getContent(), tempFile, StandardCopyOption.REPLACE_EXISTING)
                
                // Extract the zip file
                extractZip(tempFile.toString(), destinationFolder)
                
            } else {
                println "Failed to download the zip file. HTTP status code: ${response.getStatusLine().getStatusCode()}"
            }
        } finally {
            // Close HTTP client
            httpClient.close()
        }
    }

    public static void extractZip(String zipFilePath, String destinationFolder) {
        try {
            // Open zip input stream
            def zipInputStream = new ZipInputStream(new FileInputStream(zipFilePath))
            
            // Create destination folder if it does not exist
            Files.createDirectories(Paths.get(destinationFolder))
            
            // Extract each entry in the zip file
            def entry = zipInputStream.getNextEntry()
            while (entry != null) {
                def entryName = entry.getName()
                def outputFile = Paths.get(destinationFolder, entryName)
                
                // If entry is a directory, create it
                if (entry.isDirectory()) {
                    Files.createDirectories(outputFile)
                } else {
                    // Create parent directory if it does not exist
                    Files.createDirectories(outputFile.getParent())
                    
                    // Extract file
                    def outputStream = new FileOutputStream(outputFile.toFile())
                    byte[] buffer = new byte[4096]
                    def len
                    while ((len = zipInputStream.read(buffer)) > 0) {
                        outputStream.write(buffer, 0, len)
                    }
                    outputStream.close()
                }
                
                entry = zipInputStream.getNextEntry()
            }
            
            zipInputStream.closeEntry()
            zipInputStream.close()
        } catch (Exception e) {
            println "Error extracting zip file: ${e.message}"
        }
    }
    
}