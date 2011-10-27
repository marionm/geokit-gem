path = File.expand_path(File.dirname(__FILE__))
$:.unshift path unless $:.include?(path)

module Geokit
  VERSION = '1.6.2'
  # These defaults are used in Geokit::Mappable.distance_to and in acts_as_mappable
  @@default_units = :miles
  @@default_formula = :sphere

  [:default_units, :default_formula].each do |sym|
    class_eval <<-EOS, __FILE__, __LINE__
      def self.#{sym}
        if defined?(#{sym.to_s.upcase})
          #{sym.to_s.upcase}
        else
          @@#{sym}
        end
      end

      def self.#{sym}=(obj)
        @@#{sym} = obj
      end
    EOS
  end

  autoload :Inflector, 'geokit/geocoders'
  autoload :Geocoders, 'geokit/geocoders'

  autoload :Bounds,   'geokit/bounds'
  autoload :GeoLoc,   'geokit/geoloc'
  autoload :LatLng,   'geokit/latlng'
  autoload :Mappable, 'geokit/mappable'
end

# make old-style module name "GeoKit" equivalent to new-style "Geokit"
GeoKit=Geokit
