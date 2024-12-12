using System;
using System.IO;

class CRC32
{
    private static readonly uint[] CrcTable = new uint[256];

    static CRC32()
    {
        const uint polynomial = 0xEDB88320; 
        for (uint i = 0; i < 256; i++)
        {
            uint crc = i;
            for (int j = 8; j > 0; j--)
            {
                if ((crc & 1) == 1)
                    crc = (crc >> 1) ^ polynomial;
                else
                    crc >>= 1;
            }
            CrcTable[i] = crc;
        }
    }

    public static uint ComputeChecksum(byte[] data)
    {
        uint crc = 0xFFFFFFFF;

        foreach (byte b in data)
        {
            byte tableIndex = (byte)((crc ^ b) & 0xFF);
            crc = (crc >> 8) ^ CrcTable[tableIndex];
        }

        return ~crc;
    }

    public static uint ComputeChecksumForFile(string filePath)
    {
        uint crc = 0xFFFFFFFF; 

        const int bufferSize = 4096; 
        byte[] buffer = new byte[bufferSize];

        using (FileStream fs = new FileStream(filePath, FileMode.Open, FileAccess.Read))
        {
            int bytesRead;
            while ((bytesRead = fs.Read(buffer, 0, buffer.Length)) > 0)
            {
                for (int i = 0; i < bytesRead; i++)
                {
                    byte tableIndex = (byte)((crc ^ buffer[i]) & 0xFF);
                    crc = (crc >> 8) ^ CrcTable[tableIndex];
                }
            }
        }

        return ~crc; 
    }

    public static string ComputeChecksumForFileHex(string filePath)
    {
        uint crc = ComputeChecksumForFile(filePath);
        return crc.ToString("X8");
    }

    static void Main()
    {
        Console.WriteLine("Enter the full path of the file to compute its CRC32 checksum:");
        string filePath = Console.ReadLine();

        // Remove surrounding quotes if present
        filePath = filePath.Trim('"');

        // Replace single backslashes with double backslashes
        filePath = filePath.Replace("\\", "\\\\");

        if (!File.Exists(filePath))
        {
            Console.WriteLine("The file does not exist. Please check the file path.");
            return;
        }

        try
        {
            uint crc32 = ComputeChecksumForFile(filePath);
            string hexCrc32 = ComputeChecksumForFileHex(filePath);

            Console.WriteLine("File: " + filePath);
            Console.WriteLine("CRC32 Checksum (Decimal): " + crc32);
            Console.WriteLine("CRC32 Checksum (Hex): " + hexCrc32);
        }
        catch (Exception ex)
        {
            Console.WriteLine("An error occurred while computing the checksum: " + ex.Message);
        }
    }
}
