% Reading EK60 raw data file sample data
% Simrad, Lars Nonboe Andersen, 12/04-02
% mod ADR to read in raw data
% modified Gordon.Keith@csiro.au 20130528


% 	struct SampleDatagram
% 	{
% 		DatagramHeader DgHeader; // "RAW0"
% 		short Channel; // Channel number
% 		short Mode; // Datatype
% 		float TransducerDepth; // [m]
% 		float Frequency; // [Hz]
% 		float TransmitPower; // [W]
% 		float PulseLength; // [s]
% 		float BandWidth; // [Hz]
% 		float SampleInterval; // [s]
% 		float SoundVelocity; // [m/s]
% 		float AbsorptionCoefficient; // [dB/m]
% 		float Heave; // [m]
% 		float Roll; // [deg]
% 		float Pitch; // [deg]
% 		float Temperature; // [C]
% 		short TrawlUpperDepthValid; // None=0, expired=1, valid=2
% 		short TrawlOpeningValid; // None=0, expired=1, valid=2
% 		float TrawlUpperDepth; // [m]
% 		float TrawlOpening; // [m]
% 		long	  Offset; // First sample
% 		long  Count; // Number of samples
% 		short Power[]; // Compressed format - See below!
% 		short Angle[]; // See below!
% 	};
% 	
% 	power = Power * 10 * log10(2) / 256
% 	Angle - the fore-and-aft (alongship) and athwartship electrical angles
% 	are output as one 16-bit word. The alongship angle is the most significant byte
% 	while the athwartship angle is the least significant byte.
% 	Angle data is expressed in 2's complement format, and the resolution
% 	is given in steps of 180/128 electrical degrees per unit.
% 	Positive numbers denotes the fore and starboard directions

function sampledata = readsampledata(fid, detail)
sampledata.channel                  = fread(fid,1,'int16');
sampledata.mode                     = fread(fid,1,'int16');
sampledata.transducerdepth          = fread(fid,1,'float32');
sampledata.frequency                = fread(fid,1,'float32');
sampledata.transmitpower            = fread(fid,1,'float32');
sampledata.pulselength              = fread(fid,1,'float32');
sampledata.bandwidth                = fread(fid,1,'float32');
sampledata.sampleinterval           = fread(fid,1,'float32');
sampledata.soundvelocity            = fread(fid,1,'float32');
sampledata.absorptioncoefficient    = fread(fid,1,'float32');
sampledata.heave                    = fread(fid,1,'float32');
sampledata.roll                     = fread(fid,1,'float32');
sampledata.pitch                    = fread(fid,1,'float32');
sampledata.temperature              = fread(fid,1,'float32');
sampledata.trawlupperdepthvalid     = fread(fid,1,'int16');
sampledata.trawlopeningvalid        = fread(fid,1,'int16');
sampledata.trawlupperdepth          = fread(fid,1,'float32');
sampledata.trawlopening             = fread(fid,1,'float32');
sampledata.offset                   = fread(fid,1,'int32');
sampledata.count                    = fread(fid,1,'int32');

if detail
    sampledata.raw_power = fread(fid,sampledata.count,'int16');
    sampledata.power = sampledata.raw_power*10*log10(2)/256;
    
    if sampledata.mode > 1
        sampledata.raw_angle = fread(fid,[2 sampledata.count],'int8');
        sampledata.athwartship = sampledata.raw_angle(1,:)' * 180 / 128;
        sampledata.alongship   = sampledata.raw_angle(2,:)' * 180 / 128;
    end
else
    if sampledata.mode > 1
        fseek(fid,sampledata.count * 4, 0);
    else
        fseek(fid,sampledata.count * 2, 0);
    end
end



