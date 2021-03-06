module Geokit
  class LatLng
    include Mappable

    attr_accessor :lat, :lng

    # Accepts latitude and longitude or instantiates an empty instance
    # if lat and lng are not provided. Converted to floats if provided
    def initialize(lat=nil, lng=nil)
      lat = lat.to_f if lat && !lat.is_a?(Numeric)
      lng = lng.to_f if lng && !lng.is_a?(Numeric)
      @lat = lat
      @lng = lng
    end

    # Latitude attribute setter; stored as a float.
    def lat=(lat)
      @lat = lat.to_f if lat
    end

    # Longitude attribute setter; stored as a float;
    def lng=(lng)
      @lng=lng.to_f if lng
    end

    # Returns the lat and lng attributes as a comma-separated string.
    def ll
      "#{lat},#{lng}"
    end

    #returns a string with comma-separated lat,lng values
    def to_s
      ll
    end

    #returns a two-element array
    def to_a
      [lat,lng]
    end
    # Returns true if the candidate object is logically equal.  Logical equivalence
    # is true if the lat and lng attributes are the same for both objects.
    def ==(other)
      other.is_a?(LatLng) ? self.lat == other.lat && self.lng == other.lng : false
    end

    def hash
      lat.hash + lng.hash
    end

    def eql?(other)
      self == other
    end

    # A *class* method to take anything which can be inferred as a point and generate
    # a LatLng from it. You should use this anything you're not sure what the input is,
    # and want to deal with it as a LatLng if at all possible. Can take:
    #  1) two arguments (lat,lng)
    #  2) a string in the format "37.1234,-129.1234" or "37.1234 -129.1234"
    #  3) a string which can be geocoded on the fly
    #  4) an array in the format [37.1234,-129.1234]
    #  5) a LatLng or GeoLoc (which is just passed through as-is)
    #  6) anything which acts_as_mappable -- a LatLng will be extracted from it
    def self.normalize(thing,other=nil)
      # if an 'other' thing is supplied, normalize the input by creating an array of two elements
      thing=[thing,other] if other

      if thing.is_a?(String)
        thing.strip!
        if match=thing.match(/(\-?\d+\.?\d*)[, ] ?(\-?\d+\.?\d*)$/)
          return Geokit::LatLng.new(match[1],match[2])
        else
          res = Geokit::Geocoders::MultiGeocoder.geocode(thing)
          return res if res.success?
          raise Geokit::Geocoders::GeocodeError
        end
      elsif thing.is_a?(Array) && thing.size==2
        return Geokit::LatLng.new(thing[0],thing[1])
      elsif thing.is_a?(LatLng) # will also be true for GeoLocs
        return thing
      elsif thing.class.respond_to?(:acts_as_mappable) && thing.class.respond_to?(:distance_column_name)
        return thing.to_lat_lng
      elsif thing.respond_to? :to_lat_lng
        return thing.to_lat_lng
      end

      raise ArgumentError.new("#{thing} (#{thing.class}) cannot be normalized to a LatLng. We tried interpreting it as an array, string, Mappable, etc., but no dice.")
    end

    # Reverse geocodes a LatLng object using the MultiGeocoder (default), or optionally
    # using a geocoder of your choosing. Returns a new Geokit::GeoLoc object
    #
    # ==== Options
    # * :using  - Specifies the geocoder to use for reverse geocoding. Defaults to
    #             MultiGeocoder. Can be either the geocoder class (or any class that
    #             implements do_reverse_geocode for that matter), or the name of
    #             the class without the "Geocoder" part (e.g. :google)
    #
    # ==== Examples
    # LatLng.new(51.4578329, 7.0166848).reverse_geocode # => #<Geokit::GeoLoc:0x12dac20 @state...>
    # LatLng.new(51.4578329, 7.0166848).reverse_geocode(:using => :google) # => #<Geokit::GeoLoc:0x12dac20 @state...>
    # LatLng.new(51.4578329, 7.0166848).reverse_geocode(:using => Geokit::Geocoders::GoogleGeocoder) # => #<Geokit::GeoLoc:0x12dac20 @state...>
    def reverse_geocode(options = { :using => Geokit::Geocoders::MultiGeocoder })
      if options[:using].is_a?(String) or options[:using].is_a?(Symbol)
        provider = Geokit::Geocoders.const_get("#{Geokit::Inflector::camelize(options[:using].to_s)}Geocoder")
      elsif options[:using].respond_to?(:do_reverse_geocode)
        provider = options[:using]
      else
        raise ArgumentError.new("#{options[:using]} is not a valid geocoder.")
      end

      provider.send(:reverse_geocode, self)
    end
  end
end
